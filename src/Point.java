public interface Point<P extends  Point<P,V>, V extends Vector<V>> {
    P copy();
    double[] getCoordinates();
    double distanceTo(P point);
    double distanceTo(Line<P,V> line);
    V vectorTo(P point);
    int dim();
}
