package nj;

/*=======================================================================*/
/* Class: ClosestPair                                                    */
/*=======================================================================*/
/* Traverses a matrix and finds and stores the pair of taxa with the     */
/* smallest distance between them.                                       */
/*=======================================================================*/
public class ClosestPair
{
    private int f;
    private int g;
    private double qScore;

    public ClosestPair(double[][] Q)
    {
        double min = Double.POSITIVE_INFINITY;

        // Only traverse the diagonal.
        int diagonalLimit = 1;
        for (int i = 0; i < Q.length; i++)
        {
            for (int j = 0; j < diagonalLimit; j++)
            {
                if ((i != j) && (Q[i][j] < min))
                {
                    min = Q[i][j];
                    this.f = j;
                    this.g = i;
                }
            }
            diagonalLimit++;
        }

        this.qScore = Q[f][g];

    }   // End constructor ClosestPair.


    public int getF()
    {
        return this.f;
    }

    public int getG()
    {
        return this.g;
    }

    public double getQScore()
    {
        return this.qScore;
    }

    @Override
    public String toString()
    {
        return ("    Min distance pair: " +
                "f=" + this.f
                + ", g=" + this.g
                + ", q score=" + this.qScore);
    }
}
