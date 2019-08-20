package nj;

import exception.InconsistentAlignmentWidthException;
import tree.TaxonomicTreeBuilder;
import tree.TreeNode;
import util.Alignment;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;


/*=======================================================================*/
/* Class: NeighborJoining                                                */
/*=======================================================================*/
/* Implementation of the neighbor-joining algorithm based on class       */
/* lectures. The resulting trees are represented in Newick string        */
/* format.                                                               */
/*                                                                       */
/* Some considerations:                                                  */
/* - Because we are using random genera, it is possible that the         */
/*       p-distance between two sequences wlil be so great that the      */
/*       resulting Jukes-Cantor distance will attempt to take the        */
/*       natural log of a negative number, resulting in a NaN value.     */
/*       If the Jukes-Cantor distance is NaN, the program instead takes  */
/*       the p-distance as a proxy of the distance between sequences i   */
/*       and j. Otherwise, it takes the Jukes-Cantor distance.           */
/* - Because NaN Jukes-Cantor distances are replaced with the p-distance */
/*       between two sequences, it is possible that δ(f, u) or δ(g, u)   */
/*       will be a negative number. To correct this and make all branch  */
/*       lengths non-negative, the program keeps track of the negative   */
/*       branch length with the largest magnitude. When the NJ algorithm */
/*       ends, the program corrects each branch length by adding the     */
/*       magnitude of this length to all branch lengths. Therefore, the  */
/*       lowest possible branch length produced by this program is 0.    */
/*=======================================================================*/
public class NeighborJoining
{
    private static String projectDir
             = System.getProperty("user.dir");

    private static DecimalFormat df
             = new DecimalFormat("#.#####");
    
    private TaxonomicTreeBuilder taxoTree;
    private int                  numGenera;
    private boolean              sampleMode;
    private double               largestMagnitudeNegativeNumber;
    private ArrayList<Double>    branchDistances;


    /*===================================================================*/
    /* Constructor: NeighborJoining                                      */
    /*-------------------------------------------------------------------*/
    /* Initializes the parameters set by the user, the negative number   */
    /* with the largest magnitude, and the structure holding the branch  */
    /* distances.                                                        */
    /*                                                                   */
    /* Also creates a TaxonomicTreeBuilder from the .fa file built by    */
    /* CO-ARBitrator. Each node of this tree represents a record in the  */
    /* .fa file, therefore holding the information of that organism's    */
    /* taxonomic groups and accession number. This accession number is   */
    /* used to fetch the corresponding genetic sequence. Random genera   */
    /* and their sequences will be selected from this tree to use in     */
    /* each run of this program.                                         */
    /*===================================================================*/
    public NeighborJoining(int numGenera, boolean sampleMode)
    {
        this.numGenera   = numGenera;
        this.sampleMode  = sampleMode;

        this.largestMagnitudeNegativeNumber = Double.POSITIVE_INFINITY;
        this.branchDistances = new ArrayList<>();

        // Make the tree: TaxonomicTreeBuilder.
        File f = new File(projectDir + "/data/AdverbCOI.fa");
        this.taxoTree = new TaxonomicTreeBuilder(f, false, 3, true);

    }   // End constructor NeighborJoining.


    /*===================================================================*/
    /* Method: run()                                                     */
    /*-------------------------------------------------------------------*/
    /* Gets random genera from the taxonomic tree. Because the D and Q   */
    /* matrices will be reduced at every iteration of this algorithm,    */
    /* the changing indices of the records (accession numbers) are       */
    /* tracked in matrixIndexToAccessionNumMap.                          */
    /*                                                                   */
    /* Afterwards, all records are written to a .fa "training" file and  */
    /* each pair of sequences is aligned using Clustal Omega. The        */
    /* initial distance matrix (D) is created from these alignments,     */
    /* and it represents the pairwise distance between each pair of      */
    /* sequences in the training file.                                   */
    /*                                                                   */
    /* Once the initial D is created, the neighbor-joining algorithm     */
    /* runs through a loop until all taxa have been joined. At each      */
    /* iteration through the loop:                                       */
    /*     1. Create a secondary matrix, Q. The indices of Q correspond  */
    /*            to the taxa represented by the same indices in D.      */
    /*     2. From Q, a closest pair, f and g, is chosen. This pair has  */
    /*            the smallest distance between them. This pair is       *
    /*             "extruded" and their new common ancestor node is u.   */
    /*     3. The program computes the new distances from each taxon to  */
    /*            the new node, u.                                       */
    /*     4. The program creates a new D matrix, where f and g are      */
    /*        replaced with u. Therefore, the new D matrix is of size    */
    /*        (n-1) by (n-1). Using this new D matrix, the program       */
    /*        goes back to step 1.                                       */
    /*===================================================================*/
    public void run()
    {
        Map<String, String> newickSubTreeMap = new LinkedHashMap<>();

        // Choose random genera.
        List<TreeNode> genera
                = this.taxoTree.chooseRandomGenera(numGenera);

        this.printChosenGenera(genera);

        // Get sequences for all 3 genera into one map.
        // Map keys correspond to the accession number's index
        // in the D and Q matrices.
        Map<Integer, String> matrixIndexToAccessionNumMap
                = this.createMapOfMatrixIndexToAccessionNum(genera);

        String baseFileName = this.buildGeneraStringFileName(genera);

        // Build training file.
        StringBuilder trainPath = new StringBuilder(projectDir
                + "/train_files/");

        this.createDirectoryIfDoesNotExist(new File(trainPath.toString()));
        trainPath.append(baseFileName).append("_train").append(".fa");
        this.writeSequencesToTrainingFile(matrixIndexToAccessionNumMap,
                trainPath.toString());

        // Build MSA file.
        StringBuilder msaPath = new StringBuilder(projectDir
                + "/msa_files/");

        this.createDirectoryIfDoesNotExist(new File(msaPath.toString()));
        msaPath.append(baseFileName).append("_aligned.fa");
        this.writeMSAToFile(msaPath.toString(), trainPath.toString());

        try
        {
            // Create distance matrix.
            double[][] D;
            if (this.sampleMode)
            {
                matrixIndexToAccessionNumMap = this.getSampleMap();
                D = this.getSampleDistanceMatrix();
            }
            else
            {
                Alignment alignment = new Alignment(
                        new File(msaPath.toString()));

                D = this.createInitialDistanceMatrix(alignment);
            }

            // Execute until all nodes are joined.
            while (D.length > 1)
            {
                double[][] Q = this.createSecondaryMatrix(D);

                // Find closest pair. This pair is "joined"
                // to a new node, u.
                ClosestPair closestPair = new ClosestPair(Q);
                int f = closestPair.getF();
                int g = closestPair.getG();

                // Calculate new distance from closestPair to u.
                double fToU = Double.valueOf(df.format(
                        this.calculateFToUDistance(D, f, g)));

                double gToU = Double.valueOf(df.format(
                        this.calculateGToUDistance(D[f][g], fToU)));

                newickSubTreeMap
                        = this.buildNewickSubTreesAndAddBranchDistances(
                        newickSubTreeMap,
                        matrixIndexToAccessionNumMap, D,
                        f, g, fToU, gToU);

                // Calculate new distances from everything else to u.
                D = this.calculateNewDistanceTableAndDropFAndG(D, f, g);

                // After dropping F and G, re-order remaining nodes.
                matrixIndexToAccessionNumMap
                        = this.reorderRowAndColumnLabels(
                                matrixIndexToAccessionNumMap, f, g);
            }

            // Final tree is the single value of the map.
            String newickTreeString
                    = new ArrayList<>(newickSubTreeMap.values()).get(0);

            // Correct tree to have no negative branch lengths.
            if (this.largestMagnitudeNegativeNumber
                != Double.POSITIVE_INFINITY)
            {
                newickTreeString
                        = this.correctNewickTree(newickTreeString);
            }
        }
        catch (IOException | InconsistentAlignmentWidthException e)
        {
            e.printStackTrace();
        }

    }   // End run().


    /*===================================================================*/
    /* Method: createMapOfMatrixIndexToAccessionNum()                    */
    /*-------------------------------------------------------------------*/
    /* Creates initial map of matrix index to accession num. The index   */
    /* given starts at 0 and increases to match the indices in the       */
    /* matrix.                                                           */
    /*===================================================================*/
    private Map<Integer, String> createMapOfMatrixIndexToAccessionNum(
            List<TreeNode> genera)
    {
        Map<Integer, String> matrixIndexToAccessionNumMap
                = new LinkedHashMap<>();

        int index = 0;
        for (TreeNode genus : genera)
        {
            for (TreeNode species : genus.getChildren().values())
            {
                matrixIndexToAccessionNumMap.put(index,
                        species.getAccessionNum());
                index++;
            }
        }

        return matrixIndexToAccessionNumMap;

    }   // End createMapOfMatrixIndexToAccessionNum().


    /*===================================================================*/
    /* Method: buildGeneraStringFileName()                               */
    /*-------------------------------------------------------------------*/
    /* Builds base filename from chosen genera for training set and MSA  */
    /* filenames. Base filename will be in the form of:                  */
    /* "<genus1>_<genus2>_..._<genusN>".                                 */
    /*===================================================================*/
    private String buildGeneraStringFileName(List<TreeNode> genera)
    {
        StringBuilder generaSb = new StringBuilder();
        for (int i = 0; i < genera.size(); i++)
        {
            TreeNode genus = genera.get(i);

            generaSb.append(genus.getGroupName());

            if (i != genera.size() - 1)
            {
                generaSb.append("_");
            }
        }

        return generaSb.toString();

    }   // End buildGeneraStringFileName().

    
    /*===================================================================*/
    /* Method: writeSequencesToTrainingFile()                            */
    /*-------------------------------------------------------------------*/
    /* Writes all sequences from every genus into one file.              */
    /*===================================================================*/
    private void writeSequencesToTrainingFile(
            Map<Integer, String> matrixIndexToAccessionNumMap,
            String trainFileName)
    {
        try (PrintWriter trainPw = new PrintWriter(
                new BufferedWriter(new FileWriter(trainFileName))))
        {
            for (String accessionNum
                 : matrixIndexToAccessionNumMap.values())
            {
                String sequence = this.taxoTree
                        .findSequenceInFileDatabase(accessionNum);

                // Write train sequences to file.
                trainPw.print(">" + accessionNum + "\r\n");

                trainPw.print(sequence + "\r\n");
            }
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }

    }   // End writeSequencesToTrainingFile().


    /*===================================================================*/
    /* Method: writeMSAToFile()                                          */
    /*-------------------------------------------------------------------*/
    /* Externally calls Clustal Omega and writes the resulting MSA to    */
    /* file.                                                             */
    /*===================================================================*/
    private void writeMSAToFile(String MSAFileName, String trainFileName)
    {
        // Clustal doesn't overwrite old MSA files by itself.
        this.deleteFileIfExists(MSAFileName);

        try
        {
            Process aligner = new ProcessBuilder("clustalo",
                    "--infile=" + trainFileName,
                    "--outfile=" + MSAFileName)
                    .start();
            aligner.waitFor();
        }
        catch (IOException | InterruptedException e)
        {
            e.printStackTrace();
        }

    }   // End writeMSAToFile().


    /*===================================================================*/
    /* Method: getSampleMap()                                            */
    /*-------------------------------------------------------------------*/
    /* Returns initial map of matrix indices to letter labels from       */
    /* Wikipedia.                                                        */
    /* Ref: https://en.wikipedia.org/wiki/Neighbor_joining               */
    /*===================================================================*/
    private Map<Integer, String> getSampleMap()
    {
        Map<Integer, String> sampleMap = new LinkedHashMap<>();
        String alpha = "abcde";

        for (int i = 0; i < alpha.length(); i++)
        {
            sampleMap.put(i, String.valueOf(alpha.charAt(i)));
        }

        return sampleMap;

    }   // End getSampleMap().


    /*===================================================================*/
    /* Method: getSampleDistanceMatrix()                                 */
    /*-------------------------------------------------------------------*/
    /* Returns a sample initial D matrix from Wikipedia.                 */
    /* Ref: https://en.wikipedia.org/wiki/Neighbor_joining               */
    /*===================================================================*/
    private double[][] getSampleDistanceMatrix()
    {
        return new double[][] {
                {0, 5, 9, 9, 8},
                {5, 0, 10, 10, 9},
                {9, 10, 0, 8, 7},
                {9, 10, 8, 0, 3},
                {8, 9, 7, 3, 0}
        };

    }   // End getSampleDistanceMatrix().


    /*===================================================================*/
    /* Method: createInitialDistanceMatrix()                             */
    /*-------------------------------------------------------------------*/
    /* Builds the initial distance table based on the distances          */
    /* calculated using the multiple sequence alignment.                 */
    /*===================================================================*/
    private double [][] createInitialDistanceMatrix(Alignment alignment)
    {
        double [][] D = new double[alignment.size()]
                [alignment.size()];

        for (int i = 0; i < alignment.size()-1; i++)
        {
            for (int j = i+1; j < alignment.size(); j++)
            {
                if (i != j)
                {
                    String s1 = alignment.get(i);
                    String s2 = alignment.get(j);;

                    int numMismatches = this
                        .getNumMismatchesIgnoringIndelOnlyCols(s1, s2);
                    int alignmentLength = this
                        .getAlignmentLengthIgnoringIndelOnlyCols(s1, s2);
                    double p = (double) numMismatches / alignmentLength;

                    // Jukes-Cantor distance.
                    double mu = -0.75 * Math.log(1 - (4*p)/3);

                    mu = (Double.isNaN(mu) || Double.isInfinite(mu))
                          ? p : mu;
                    mu = Double.valueOf(df.format(mu));

                    // Distance matrix is symmetric.
                    D[i][j] = mu;
                    D[j][i] = mu;
                }
            }
        }

        return D;

    }   // End createInitialDistanceMatrix().


    /*===================================================================*/
    /* Method: getNumMismatchesIgnoringIndelOnlyCols()                   */
    /*-------------------------------------------------------------------*/
    /* Gets total number of mismatches between two sequences from a      */
    /* multiple sequence alignment, ignoring indel-only columns. This    */
    /* number is used to calculate p-distance between two sequences.     */
    /*===================================================================*/
    private int getNumMismatchesIgnoringIndelOnlyCols(
            String s1, String s2)
    {
        int numMismatchesIgnoringIndelOnlyCols = 0;
        for (int i = 0; i < s1.length(); i++)
        {
            if (!((s1.charAt(i) == '-') && (s2.charAt(i) == '-')) &&
                    (s1.charAt(i) != s2.charAt(i)))
            {
                numMismatchesIgnoringIndelOnlyCols++;
            }
        }

        return numMismatchesIgnoringIndelOnlyCols;

    }   // End getNumMismatchesIgnoringIndelOnlyCols().


    /*===================================================================*/
    /* Method: getAlignmentLengthIgnoringIndelOnlyCols()                 */
    /*-------------------------------------------------------------------*/
    /* Gets total length of the alignment, ignoring indel-only columns.  */
    /* Used to calculate p-distance between two sequences.               */
    /*===================================================================*/
    private int getAlignmentLengthIgnoringIndelOnlyCols(
            String s1, String s2)
    {
        int length = s1.length();
        for (int i = 0; i < s1.length(); i++)
        {
            if ((s1.charAt(i) == '-') && (s2.charAt(i) == '-'))
            {
                length--;
            }
        }

        return length;

    }   // End getAlignmentLengthIgnoringIndelOnlyCols().


    /*===================================================================*/
    /* Method: createSecondaryMatrix()                                   */
    /*-------------------------------------------------------------------*/
    /* Builds the neighbor-joining Q matrix.                             */
    /*===================================================================*/
    private double [][] createSecondaryMatrix(double [][] D)
    {
        // Same number and columns as D.
        double [][] Q = new double[D.length][D[0].length];

        int n = D.length;
        int diagonalLimit = 1;

        for (int i = 0; i < Q.length; i++)
        {
            for (int j = 0; j < diagonalLimit; j++)
            {
                if (i != j)
                {
                    // Given equations.
                    double distanceFromIToEverything
                            = this.calculateDistanceToEverything(D, i);

                    double distanceFromJToEverything
                            = this.calculateDistanceToEverything(D, j);

                    double secondaryDistance
                            = Double.valueOf(df.format(
                                    (n - 2) * D[i][j]
                                            - distanceFromIToEverything
                                            - distanceFromJToEverything));

                    // Set new distance.
                    Q[i][j] = secondaryDistance;
                    Q[j][i] = secondaryDistance;
                }
                else
                {
                    Q[i][j] = 0;
                    Q[j][i] = 0;
                }

            }
            diagonalLimit++;
        }

        return Q;

    }   // End createSecondaryMatrix().


    /*===================================================================*/
    /* Method: calculateDistanceToEverything()                           */
    /*-------------------------------------------------------------------*/
    /* Calculates the distance from the current taxon to all the other   */
    /* taxa in the distance matrix.                                      */
    /*===================================================================*/
    private double calculateDistanceToEverything(double [][] D,
                                                 int current)
    {
        double distance = 0;

        for (int j = 0; j < D.length; j++)
        {
            if (j != current)
            {
                distance += D[current][j];
            }
        }

        return distance;

    }   // End calculateDistanceToEverything().


    /*===================================================================*/
    /* Method: calculateFToUDistance()                                   */
    /*-------------------------------------------------------------------*/
    /* Given equation.                                                   */
    /*===================================================================*/
    private double calculateFToUDistance(double[][] D, int f, int g)
    {
        return (0.5 * D[f][g]) + (1.0 / (2 * (D.length-2)))
                * (this.calculateDistanceToEverything(D, f)
                - this.calculateDistanceToEverything(D, g));

    }   // End calculateFToUDistance().


    /*===================================================================*/
    /* Method: calculateGToUDistance()                                   */
    /*-------------------------------------------------------------------*/
    /* Given equation.                                                   */
    /*===================================================================*/
    private double calculateGToUDistance(double fToG, double fToU)
    {
        return fToG - fToU;

    }   // End calculateGToUDistance().

    
    /*===================================================================*/
    /* Method: buildNewickSubTreesAndAddBranchDistances()                */
    /*-------------------------------------------------------------------*/
    /* Stores the Newick trees for all joined taxa. Taxa can be joined   */
    /* at different rates; for example, A and B can be joined to A_B,    */
    /* then C and D can be joined to C_D, before A_B and C_D are joined. */
    /* Therefore, the program needs to store the separate A_B and C_D    */
    /* Newick trees. Once A_B and C_D are joined, it can combine their   */
    /* Newick trees to A_B_C_D.                                          */
    /*                                                                   */
    /* This method keeps track of joined node names as the key and their */
    /* corresponding Newick trees as the values. Whenever the method is  */
    /* called, it creates a new subtree map. If f and g already exist    */
    /* as keys in the previous map, this method joins their Newick       */
    /* trees. If only f exists and g does not, it joins the Newick tree  */
    /* of f with the g node. If both do not exist in the previous map,   */
    /* this method creates a new key with a new Newick tree:             */
    /*     f_g : (f: fToU, g: gToU).                                     */
    /*                                                                   */
    /* In any case, a new key is created. If f and/or g already exist in */
    /* the previous map, this method excludes them from the new map      */
    /* because they have become part of a new node.                      */
    /*                                                                   */
    /* Finally, this method copies over all unaffected elements from the */
    /* previous map. In the above example, if the previous map contained */
    /* the key E_F, this key and value would be directly copied to the   */
    /* new map because it was not involved in the joining of the new     */
    /* node.                                                             */
    /*===================================================================*/
    private Map<String, String> buildNewickSubTreesAndAddBranchDistances(
            Map<String, String> prevMap,
            Map<Integer, String> matrixIndexToAccessionNumMap,
            double[][] D, int f, int g, double fToU, double gToU)
    {
        Map<String, String> newSubTreeMap = new LinkedHashMap<>();

        String fName = matrixIndexToAccessionNumMap.get(f);
        String gName = matrixIndexToAccessionNumMap.get(g);
        String combinedName = fName + "__" + gName;

        StringBuilder newSubTreeSb = new StringBuilder();

        if (!Double.isNaN(fToU))
        {
            if (!fName.contains("__"))
            {
                newSubTreeSb.append("(").append(fName)
                        .append(": ").append(fToU);
            }
            else if (fName.contains("__"))
            {
                // fName already has a subtree.
                newSubTreeSb.append("(").append(prevMap.get(fName))
                        .append(": ").append(fToU);
            }
            this.branchDistances.add(fToU);
            this.checkIfScoreIsLargestNegativeNum(fToU);
        }
        else
        {
            if (fName.contains("__"))
            {
                newSubTreeSb.append(prevMap.get(fName))
                        .setLength(newSubTreeSb.length()-1);
            }
        }

        double gDistance = gToU;
        if (Double.isNaN(gToU))
        {
            gDistance = D[f][g];
        }
        if (!gName.contains("__"))
        {
            newSubTreeSb.append(", ").append(gName)
                    .append(": ").append(gDistance).append(")");
        }
        else if (gName.contains("__"))
        {
            // gName already has a subtree.
            newSubTreeSb.append(", ").append(prevMap.get(gName))
                    .append(": ").append(gDistance).append(")");
        }
        this.branchDistances.add(gDistance);
        this.checkIfScoreIsLargestNegativeNum(gDistance);

        // Add new subtree.
        newSubTreeMap.put(combinedName, newSubTreeSb.toString());

        // Add all previous subtrees excluding f and g.
        for (String key : prevMap.keySet())
        {
            if ((key.compareToIgnoreCase(fName) != 0) &&
                    (key.compareToIgnoreCase(gName) != 0))
            {
                newSubTreeMap.put(key, prevMap.get(key));
            }
        }
        return newSubTreeMap;

    }   // End buildNewickSubTreesAndAddBranchDistances().


    /*===================================================================*/
    /* Method: checkIfScoreIsLargestNegativeNum()                        */
    /*-------------------------------------------------------------------*/
    /* If score is negative and the magnitude of the score is greater    */
    /* than the current largest magnitude negative number, the score     */
    /* becomes the new largest magnitude negative number. This is used   */
    /* to correct the Newick tree so that there are no negative branch   */
    /* lengths.                                                          */
    /*===================================================================*/
    private void checkIfScoreIsLargestNegativeNum(double score)
    {

        if ((score < 0) && (score < this.largestMagnitudeNegativeNumber))
        {
            this.largestMagnitudeNegativeNumber = score;
        }

    }   // End checkIfScoreIsLargestNegativeNum().


    /*===================================================================*/
    /* Method: calculateNewDistanceTableAndDropFAndG()                   */
    /*-------------------------------------------------------------------*/
    /* Calculates new distances to the common ancestor of f and g.       */
    /* The new common ancestor node (u) has the row and column index 0.  */
    /* All other taxa besides f and g are shifted down so that the new   */
    /* distance table excludes f and g.                                  */
    /*===================================================================*/
    private double[][] calculateNewDistanceTableAndDropFAndG(
            double[][] oldD, int f, int g)
    {
        double[][] newD = new double[oldD.length-1][oldD[0].length-1];

        int newDRow = 1;
        int newDCol = 1;

        int uIndex = 0;
        double[] newUDistances = new double[newD.length];
        newUDistances[uIndex++] = 0.0;

        for (int i = 0; i < oldD.length; i++)
        {
            if ((i == f) || (i == g))
            {
                continue;
            }
            for (int j = 0; j < oldD.length; j++)
            {
                // Copy over all indices, except for
                // rows/columns f and g, from oldD to newD.
                if ((j == f) || (j == g))
                {
                    continue;
                }
                newD[newDRow][newDCol] = oldD[i][j];
                newDCol++;

            }
            newUDistances[uIndex++] = Double.valueOf(
                    df.format(0.5 * (oldD[f][i] + oldD[g][i]
                    - oldD[f][g])));
            newDRow++;
            newDCol = 1;
        }

        for (int i = 0; i < newD.length; i++)
        {
            newD[i][0] = newUDistances[i];
            newD[0][i] = newUDistances[i];
        }

        return newD;

    }   // End calculateNewDistanceTableAndDropFAndG().


    /*===================================================================*/
    /* Method: reorderRowAndColumnLabels()                               */
    /*-------------------------------------------------------------------*/
    /* Re-assigns keys to values to update the location of each          */
    /* accession number in the matrices. The newly joined pair are given */
    /* index 0. The remaining values from the old map are added in the   */
    /* same order following the newly joined pair.                       */
    /*===================================================================*/
    private Map<Integer, String> reorderRowAndColumnLabels(
            Map<Integer, String> oldMap, int f, int g)
    {
        Map<Integer, String> newMap = new LinkedHashMap<>();

        // Index 0 is the new joined group.
        int newIndex = 0;

        // Represent joined taxa as "taxon1__taxon2".
        newMap.put(newIndex++, oldMap.get(f) + "__" + oldMap.get(g));

        for (Integer index : oldMap.keySet())
        {
            if ((index != f) && (index != g))
            {
                newMap.put(newIndex++, oldMap.get(index));
            }
        }

        return newMap;

    }   // End reorderRowAndColumnLabels().


    /*===================================================================*/
    /* Method: correctNewickTree()                                       */
    /*-------------------------------------------------------------------*/
    /* Corrects the Newick tree string by adding the magnitude of the    */
    /* smallest negative number to all the original branch lengths, then */
    /* replacing the branch lengths in the string with the new, all non- */
    /* negative, branch lengths.                                         */
    /*===================================================================*/
    private String correctNewickTree(String newickTree)
    {
        String correctedTree = newickTree;
        double correction = Math.abs(this.largestMagnitudeNegativeNumber);

        for (Double oldDistance : this.branchDistances)
        {
            Double newDistance = oldDistance + correction;

            DecimalFormat df = new DecimalFormat("#.#####");
            newDistance = Double.valueOf(df.format(newDistance));

            correctedTree = correctedTree.replaceAll(
                    String.valueOf(oldDistance),
                    String.valueOf(newDistance));
        }

        return correctedTree;

    }   // End correctNewickTree().


    /*===================================================================*/
    /* Method: printChosenGenera()                                       */
    /*-------------------------------------------------------------------*/
    /* Prints each random genus that was chosen, along with its size     */
    /* (the number of species in that genus).                            */
    /*===================================================================*/
    private void printChosenGenera(List<TreeNode> genera)
    {
        for (int i = 0; i < genera.size(); i++)
        {
            TreeNode genus = genera.get(i);
        }

    }   // End printChosenGenera().


    /*===================================================================*/
    /* Method: createDirectoryIfDoesNotExist(File)                       */
    /*===================================================================*/
    private void createDirectoryIfDoesNotExist(File dir)
    {
        if (!dir.exists())
        {
            dir.mkdir();
        }

    }   // End createDirectoryIfDoesNotExist().


    /*===================================================================*/
    /* Method: createDirectoryIfDoesNotExist(String)                     */
    /*===================================================================*/
    private void createDirectoryIfDoesNotExist(String dir)
    {
        this.createDirectoryIfDoesNotExist(new File(dir));

    }   // End createDirectoryIfDoesNotExist().


    /*===================================================================*/
    /* Method: deleteFileIfExists(File)                                  */
    /*===================================================================*/
    private void deleteFileIfExists(File file)
    {
        if (file.exists()) { file.delete(); }

    }   // End deleteFileIfExists(File).


    /*===================================================================*/
    /* Method: deleteFileIfExists(String)                                */
    /*===================================================================*/
    private void deleteFileIfExists(String fileName)
    {
        this.deleteFileIfExists(new File(fileName));

    }   // End deleteFileIfExists(String).


    /*===================================================================*/
    /* Method: matrixToString()                                          */
    /*-------------------------------------------------------------------*/
    /* Returns a matrix in String format.                                */
    /*===================================================================*/
    private String matrixToString(double [][] matrix)
    {
        StringBuilder s = new StringBuilder();
        int maxLength = 0;

        int diagonalLimit = 1;
        for (int i = 0; i < matrix.length; i++)
        {
            for (int j = 0; j < diagonalLimit; j++)
            {
                if (String.valueOf(matrix[i][j]).length() > maxLength)
                {
                    maxLength = String.valueOf(matrix[i][j]).length();
                }
            }
            diagonalLimit++;
        }

        String maxLengthString = "%" + maxLength + "s";
        for (int i = 0; i < matrix.length; i++)
        {
            for (int j = 0; j < matrix[0].length; j++)
            {
                s.append(String.format(maxLengthString,
                                       String.valueOf(matrix[i][j])));

                if (j != matrix[0].length-1)
                {
                    s.append("    ");
                }
            }
            s.append("\n");
        }

        return s.toString();

    }   // End matrixToString().


    /*===================================================================*/
    /* Method: mapToString()                                             */
    /*-------------------------------------------------------------------*/
    /* Returns a matrix in String format where each line contains:       */
    /* key: value.                                                       */
    /*===================================================================*/
    private String mapToString(Map<Integer, String> map)
    {
        StringBuilder s = new StringBuilder();

        for (Integer i : map.keySet())
        {
            s.append(i).append(": ").append(map.get(i)).append("\n");
        }

        return s.toString();

    }   // End mapToString().


    /*===================================================================*/
    /* Method: newickMapToString()                                       */
    /*-------------------------------------------------------------------*/
    /* Returns a matrix in String format where each line contains:       */
    /* key: value.                                                       */
    /*===================================================================*/
    private String newickMapToString(Map<String, String> map)
    {
        StringBuilder s = new StringBuilder();

        for (String key : map.keySet())
        {
            s.append(key).append(": ").append(map.get(key)).append("\n");
        }

        return s.toString();

    }   // End newickMapToString().


    /*===================================================================*/
    /* Method: main()                                                    */
    /*===================================================================*/
    public static void main(String[] args)
    {
        NeighborJoining nj = new NeighborJoining(3, true);
        nj.run();

    }   // End main().

}   // End class NeighborJoining

