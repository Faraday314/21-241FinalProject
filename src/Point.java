public interface Point<P extends  Point<P,V>, V extends math.Vector<V>> {
    P copy();
    double[] getCoordinates();
    double distanceTo(P point);
    double distanceTo(math.Line<P,V> line);
    V vectorTo(P point);
    int dim();
}
