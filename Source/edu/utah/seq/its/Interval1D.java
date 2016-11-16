package edu.utah.seq.its;

/******************************************************************************
 *  Compilation:  javac Interval1D.java
 *  Execution:    java Interval1D
 *  
 *  Interval data type with integer coordinates.
 *
 * Downloaded from http://algs4.cs.princeton.edu/home/ Robert Sedgewick and Kevin Wayne
 ******************************************************************************/


public class Interval1D implements Comparable<Interval1D>, java.io.Serializable {
    public final int min;  // min endpoint
    public final int max;  // max endpoint
    private static final long serialVersionUID = 1L;

    // precondition: min <= max
    public Interval1D(int min, int max) {
        if (min <= max) {
            this.min = min;
            this.max = max;
        }
        else throw new RuntimeException("Illegal interval");
    }

    // does this interval intersect that one?
    public boolean intersects(Interval1D that) {
        if (that.max < this.min) return false;
        if (this.max < that.min) return false;
        return true;
    }

    // does this interval a intersect b?
    public boolean contains(int x) {
        return (min <= x) && (x <= max);
    }

    public int compareTo(Interval1D that) {
        if      (this.min < that.min) return -1;
        else if (this.min > that.min) return +1;
        else if (this.max < that.max) return -1;
        else if (this.max > that.max) return +1;
        else                          return  0;
    }

    public String toString() {
        return "[" + min + ", " + max + "]";
    }




    // test client
    public static void main(String[] args) {
        Interval1D a = new Interval1D(15, 20);
        Interval1D b = new Interval1D(25, 30);
        Interval1D c = new Interval1D(10, 40);
        Interval1D d = new Interval1D(40, 50);

        System.out.println("a = " + a);
        System.out.println("b = " + b);
        System.out.println("c = " + c);
        System.out.println("d = " + d);

        System.out.println("b intersects a = " + b.intersects(a));
        System.out.println("a intersects b = " + a.intersects(b));
        System.out.println("a intersects c = " + a.intersects(c));
        System.out.println("a intersects d = " + a.intersects(d));
        System.out.println("b intersects c = " + b.intersects(c));
        System.out.println("b intersects d = " + b.intersects(d));
        System.out.println("c intersects d = " + c.intersects(d));

    }

}


