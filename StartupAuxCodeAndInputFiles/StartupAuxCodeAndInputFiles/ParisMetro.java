// Paris (P2- Part B)
// Startup code given in the Fall 2025 for csi2110/csi2510
// Completed by: Rayyan Lodhi
//
// This program reads the Paris Metro description from standard input,
// builds the directed metro graph (metroGraph), compresses it into a simpler
// undirected graph where each vertex corresponds to a "Hub Station"
// (expensiveParisSubway), and then applies Kruskal's algorithm using
// the Partition ADT to find a minimum spanning tree of the Hub Stations.
//
// The cost of each edge in expensiveParisSubway is the minimum travel time
// (in seconds) along the metro lines between two Hub Stations. Once the MST
// is computed, the program outputs the total yearly cost Peter has to pay
// for hub-to-hub passes and lists the segments (Hub Station pairs) he must buy.
//

import java.util.*;
import java.io.*;

public class ParisMetro {

    /**
     * Program entry point.
     * Reads the metro description and prints the required MST information.
     */
    public static void main(String[] args) {
        readMetro();
    }

    /**
     * Reads the metro description, constructs metroGraph and expensiveParisSubway,
     * and runs Kruskal to obtain the minimum-cost set of hub segments.
     */
    public static void readMetro() {

        Scanner inputScanner = new Scanner(System.in);

        /*
         * Phase 1: Read basic graph sizes for metroGraph.
         *
         * numVertices = number of vertices (stations on particular lines)
         * numEdges = number of edges (directed connections or walking links)
         */
        int numVertices = inputScanner.nextInt();
        int numEdges = inputScanner.nextInt();

        /*
         * Arrays to store, for each vertex index:
         *
         * vertexNumber[index] : numeric station code from the input (e.g., 0016)
         * stationName[index] : station name string (e.g., "Bastille")
         */
        int[] vertexNumber = new int[numVertices];
        String[] stationName = new String[numVertices];

        /*
         * Adjacency list for the directed metroGraph.
         * Each entry metroAdjacencyList[u] stores the outgoing MetroEdge objects
         * from vertex u with positive travel time.
         */
        @SuppressWarnings("unchecked")
        ArrayList<MetroEdge>[] metroAdjacencyList = new ArrayList[numVertices];
        for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
            metroAdjacencyList[vertexIndex] = new ArrayList<MetroEdge>();
        }

        /*
         * Maps the numeric vertex code (from the file) to the internal index
         * used in arrays (0..numVertices-1).
         */
        HashMap<Integer, Integer> vertexIdToIndex = new HashMap<Integer, Integer>();

        /*
         * Partition used to group vertices that belong to the same physical
         * station via walking edges (weight -1). These groups will become
         * Hub Stations in expensiveParisSubway.
         */
        Partition<Integer> hubPartition = new Partition<Integer>();

        /*
         * For each vertex index, vertexNodes[index] stores the Node<Integer>
         * returned by hubPartition.makeCluster when the vertex is first created.
         */
        @SuppressWarnings("unchecked")
        Node<Integer>[] vertexNodes = new Node[numVertices];

        /*
         * Phase 1a: Read all vertices.
         *
         * Each vertex line consists of:
         * - an integer vertex ID
         * - a station name (may contain spaces, so we read the rest of the line)
         */
        for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
            int vertexId = inputScanner.nextInt();
            String nameOfStation = inputScanner.nextLine().trim();

            vertexNumber[vertexIndex] = vertexId;
            stationName[vertexIndex] = nameOfStation;

            vertexIdToIndex.put(vertexId, vertexIndex);

            /*
             * Initially, each vertex is in its own cluster in the Partition.
             * Walking edges (weight -1) will later merge these clusters.
             */
            vertexNodes[vertexIndex] = hubPartition.makeCluster(Integer.valueOf(vertexIndex));
        }

        /*
         * The next token after the vertex list is the symbol "$",
         * which indicates the end of the vertex section.
         */
        inputScanner.next();

        /*
         * Phase 1b: Read all edges.
         *
         * For each edge line:
         * - id1, id2 : numeric vertex codes
         * - w : weight
         *
         * If w == -1, the edge is a walking connection and is used only to
         * merge clusters in hubPartition.
         * If w > 0, the edge is a directed metro connection in metroGraph.
         */
        for (int edgeIndex = 0; edgeIndex < numEdges; edgeIndex++) {
            int idFrom = inputScanner.nextInt();
            int idTo = inputScanner.nextInt();
            int weight = inputScanner.nextInt();

            int sourceIndex = vertexIdToIndex.get(idFrom);
            int destinationIndex = vertexIdToIndex.get(idTo);

            if (weight == -1) {
                hubPartition.union(vertexNodes[sourceIndex], vertexNodes[destinationIndex]);
            } else if (weight > 0) {
                metroAdjacencyList[sourceIndex].add(new MetroEdge(destinationIndex, weight));
            }
        }

        inputScanner.close();

        System.out.println("Paris Metro Graph has " + numVertices +
                " vertices and " + numEdges + " edges.");

        /*
         * Phase 2: Identify Hub Stations.
         *
         * A Hub Station is a physical station that contains at least two
         * different vertices (stations on possibly different lines).
         *
         * We walk over all vertices, locate their cluster leader in hubPartition,
         * and group leaders that correspond to clusters of size >= 2.
         */

        HashMap<Node<Integer>, Integer> leaderToHubIndex = new HashMap<Node<Integer>, Integer>();
        ArrayList<String> hubStationNames = new ArrayList<String>();
        ArrayList<ArrayList<Integer>> hubStationMembers = new ArrayList<ArrayList<Integer>>();

        /*
         * For each metro vertex we store the index of its Hub Station, or -1 if
         * it does not belong to any Hub Station (cluster size == 1).
         */
        int[] hubIndexOfVertex = new int[numVertices];
        Arrays.fill(hubIndexOfVertex, -1);

        for (int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++) {
            Node<Integer> leader = hubPartition.find(vertexNodes[vertexIndex]);
            int clusterSize = hubPartition.clusterSize(leader);

            if (clusterSize >= 2) {
                Integer hubIndex = leaderToHubIndex.get(leader);
                if (hubIndex == null) {
                    hubIndex = leaderToHubIndex.size();
                    leaderToHubIndex.put(leader, hubIndex);

                    /*
                     * Use the first encountered station name in the cluster
                     * as the name of this Hub Station.
                     */
                    hubStationNames.add(stationName[vertexIndex]);
                    hubStationMembers.add(new ArrayList<Integer>());
                }
                hubStationMembers.get(hubIndex).add(vertexIndex);
                hubIndexOfVertex[vertexIndex] = hubIndex;
            }
        }

        int numberOfHubs = hubStationNames.size();
        int totalHubVertices = 0;
        for (ArrayList<Integer> memberList : hubStationMembers) {
            totalHubVertices += memberList.size();
        }

        System.out.print("Hub Stations = [ ");
        for (int hubIndex = 0; hubIndex < numberOfHubs; hubIndex++) {
            if (hubIndex > 0) {
                System.out.print(", ");
            }
            System.out.print(hubStationNames.get(hubIndex));
        }
        System.out.println(" ]");
        System.out.println("Number of Hub Stations = " + numberOfHubs +
                " (total Hub Vertices = " + totalHubVertices + ")");

        /*
         * Phase 3: Build expensiveParisSubway.
         *
         * For each Hub Station h, we run a multi-source Dijkstra starting from
         * all vertices belonging to h with initial distance 0. We explore the
         * directed metroGraph until we reach another Hub Station h2. The first
         * such reach defines a candidate undirected segment between h and h2
         * with cost equal to the minimum travel time in seconds.
         *
         * bestSegmentCost[h][h2] will store the minimum cost from Hub h to Hub h2
         * found in that way. Only finite values become edges in the hub graph.
         */

        final int INF = 1000000000;
        int[][] bestSegmentCost = new int[numberOfHubs][numberOfHubs];
        for (int hubIndex = 0; hubIndex < numberOfHubs; hubIndex++) {
            Arrays.fill(bestSegmentCost[hubIndex], INF);
        }

        for (int hubIndex = 0; hubIndex < numberOfHubs; hubIndex++) {

            int[] distanceFromHub = new int[numVertices];
            Arrays.fill(distanceFromHub, INF);

            PriorityQueue<NodeDist> priorityQueue = new PriorityQueue<NodeDist>(
                    Comparator.comparingInt(node -> node.dist));

            /*
             * Initialize the queue with all vertices in this Hub Station.
             * Each is a source with distance 0.
             */
            for (int vertexIndex : hubStationMembers.get(hubIndex)) {
                distanceFromHub[vertexIndex] = 0;
                priorityQueue.add(new NodeDist(vertexIndex, 0));
            }

            /*
             * Standard Dijkstra loop, but we stop exploring when we reach
             * a different Hub Station. That reach defines the cost of a
             * hub-to-hub segment.
             */
            while (!priorityQueue.isEmpty()) {
                NodeDist current = priorityQueue.poll();
                int currentVertex = current.v;
                int currentDistance = current.dist;

                if (currentDistance != distanceFromHub[currentVertex]) {
                    continue;
                }

                for (MetroEdge outgoingEdge : metroAdjacencyList[currentVertex]) {
                    if (outgoingEdge.weight <= 0) {
                        continue;
                    }

                    int neighborVertex = outgoingEdge.to;
                    int newDistance = currentDistance + outgoingEdge.weight;

                    if (newDistance >= distanceFromHub[neighborVertex]) {
                        continue;
                    }

                    int neighborHubIndex = hubIndexOfVertex[neighborVertex];

                    /*
                     * If the neighbor is not in any Hub or is still within the same
                     * Hub, we keep exploring. Otherwise, we have reached a new Hub.
                     */
                    if (neighborHubIndex == -1 || neighborHubIndex == hubIndex) {
                        distanceFromHub[neighborVertex] = newDistance;
                        priorityQueue.add(new NodeDist(neighborVertex, newDistance));
                    } else {
                        if (newDistance < bestSegmentCost[hubIndex][neighborHubIndex]) {
                            bestSegmentCost[hubIndex][neighborHubIndex] = newDistance;
                        }
                    }
                }
            }
        }

        /*
         * Phase 4: Collect undirected segments between Hub Stations.
         *
         * For each unordered pair {i, j}, if we have a finite cost either from
         * i to j or from j to i, we take the minimum of the two directions as
         * the cost of the undirected segment.
         *
         * Each segment is stored as a small int[] of the form {i, j, weight}.
         */
        ArrayList<int[]> hubSegments = new ArrayList<int[]>();
        for (int hubIndex1 = 0; hubIndex1 < numberOfHubs; hubIndex1++) {
            for (int hubIndex2 = hubIndex1 + 1; hubIndex2 < numberOfHubs; hubIndex2++) {
                int costForward = bestSegmentCost[hubIndex1][hubIndex2];
                int costBackward = bestSegmentCost[hubIndex2][hubIndex1];
                int segmentCost = Math.min(costForward, costBackward);

                if (segmentCost < INF) {
                    hubSegments.add(new int[] { hubIndex1, hubIndex2, segmentCost });
                }
            }
        }

        System.out.println("Number of Possible Segments = " + hubSegments.size());

        /*
         * Phase 5: Run Kruskal's algorithm on expensiveParisSubway using
         * the Partition ADT from Part A.
         */

        int numberOfSegments = hubSegments.size();
        if (numberOfHubs == 0 || numberOfSegments == 0) {
            System.out.println("Impossible");
            return;
        }

        /*
         * Create one cluster per Hub Station. hubNodes[h] is the Node<String>
         * corresponding to Hub Station h in the Partition used for Kruskal.
         */
        Partition<String> hubMSTPartition = new Partition<String>();
        @SuppressWarnings("unchecked")
        Node<String>[] hubNodes = new Node[numberOfHubs];
        for (int hubIndex = 0; hubIndex < numberOfHubs; hubIndex++) {
            hubNodes[hubIndex] = hubMSTPartition.makeCluster(hubStationNames.get(hubIndex));
        }

        /*
         * Convert each hub segment into an Edge<String> object whose weight is
         * the yearly cost (in seconds) of that segment.
         */
        @SuppressWarnings("unchecked")
        Edge<String>[] mstEdgeArray = new Edge[numberOfSegments];
        for (int segmentIndex = 0; segmentIndex < numberOfSegments; segmentIndex++) {
            int[] segment = hubSegments.get(segmentIndex);
            int hubIndex1 = segment[0];
            int hubIndex2 = segment[1];
            int segmentCost = segment[2];

            mstEdgeArray[segmentIndex] = new Edge<String>(hubNodes[hubIndex1], hubNodes[hubIndex2], segmentCost);
        }

        /*
         * Sort all candidate segments by cost in non-decreasing order.
         */
        Arrays.sort(mstEdgeArray);

        /*
         * Apply Kruskal's algorithm.
         * We repeatedly add the cheapest edge that connects two different
         * components in the Partition until we have a spanning tree
         * or we run out of edges.
         */
        ArrayList<Edge<String>> mstEdges = new ArrayList<Edge<String>>();
        int totalMSTCost = 0;
        int mstEdgesUsed = 0;

        for (Edge<String> segmentEdge : mstEdgeArray) {
            Node<String> hubNode1 = segmentEdge.getU();
            Node<String> hubNode2 = segmentEdge.getV();

            Node<String> leader1 = hubMSTPartition.find(hubNode1);
            Node<String> leader2 = hubMSTPartition.find(hubNode2);

            if (leader1 != leader2) {
                hubMSTPartition.union(hubNode1, hubNode2);
                mstEdges.add(segmentEdge);
                totalMSTCost += segmentEdge.getWeight();
                mstEdgesUsed++;

                if (mstEdgesUsed == numberOfHubs - 1) {
                    break;
                }
            }
        }

        /*
         * If we could not connect all Hub Stations with numberOfHubs - 1 edges,
         * expensiveParisSubway is disconnected and the task is impossible.
         */
        if (mstEdgesUsed != numberOfHubs - 1) {
            System.out.println("Impossible");
            return;
        }

        /*
         * Phase 6: Output the MST result:
         *
         * - Total yearly cost.
         * - List of segments (Hub Station pairs) Peter must buy.
         */
        System.out.println();
        System.out.println("Total Cost = $" + totalMSTCost);
        System.out.println("Segments to Buy:");

        int segmentNumber = 1;
        for (Edge<String> mstSegment : mstEdges) {
            String hubName1 = mstSegment.getU().getElement();
            String hubName2 = mstSegment.getV().getElement();
            int cost = mstSegment.getWeight();

            System.out.println(segmentNumber + "( " + hubName1 + " - " + hubName2 + " ) - $" + cost);
            segmentNumber++;
        }
    }
}

/*
 * The following classes are the Partition ADT and Edge class given in Part A.
 * They are kept unchanged (except for formatting) so that the same data
 * structure and operations are used in both parts of the assignment.
 */

// Given in programming assignment 1.
class Cluster<E> {

    private Cluster<E> prevCluster; // pointer to previous cluster
    private Cluster<E> nextCluster; // pointer to the head of this cluster

    private Node<E> head;
    private Node<E> tail;
    /** Number of elements in the list */
    private int size = 0;

    /** Constructs a new empty list. */
    public Cluster(Cluster<E> prevCluster, Cluster<E> nextCluster) {
        head = null; // head
        tail = head; // tail
        this.prevCluster = prevCluster;
        this.nextCluster = nextCluster;
    }

    public void setNextCluster(Cluster<E> c) {
        this.nextCluster = c;
    }

    public void setPrevCluster(Cluster<E> c) {
        this.prevCluster = c;
    }

    public Cluster<E> getNextCluster() {
        return this.nextCluster;
    }

    public Cluster<E> getPrevCluster() {
        return this.prevCluster;
    }

    /**
     * Returns the number of elements in the list.
     * 
     * @return number of elements in the list
     */
    public int size() {
        return size;
    };

    /**
     * Tests whether the list is empty.
     * 
     * @return true if the list is empty, false otherwise
     */
    public boolean isEmpty() {
        return size == 0;
    };

    /**
     * Returns the first Node in the list.
     *
     * @return the first Node in the list (or null, if empty)
     */
    public Node<E> first() {
        if (this.isEmpty()) {
            return null;
        }
        return head;
    };

    /**
     * Returns the last Node in the list.
     *
     * @return the last Node in the list (or null, if empty)
     */
    public Node<E> last() {
        if (this.isEmpty()) {
            return null;
        }
        return tail;
    }

    /**
     * Returns the Node immediately after Position p.
     * 
     * @param p a Node of the list
     * @return the subsequent Node
     */

    public Node<E> after(Node<E> p) {
        return p.getNext();
    }

    /**
     * Appends a cluster c to the end of this Cluster
     * 
     * @param c Cluster to be appended
     */
    public void add(Cluster<E> c) {
        Node<E> n = c.first();
        n.setPrev(tail);
        tail.setNext(n);
        this.tail = c.tail;
        size = size + c.size;

    }

    /**
     * Inserts an element at the back of the list.
     *
     * @param e the new element
     * @return the Node representing the location of the new element
     */

    public Node<E> addLast(E e) {
        Node<E> n = new Node(null, e, null, this); // (prev, element, next, )
        if (this.isEmpty()) {
            head = n;
            tail = head;
            size++;
            return n;
        }

        n.setPrev(tail);
        tail.setNext(n);
        tail = n;
        size++;
        return n;
    }

    // Debugging code
    /**
     * Produces a string representation of the contents of the list.
     * This exists for debugging purposes only.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("(");
        Node<E> walk = head;
        while (walk != null) {
            sb.append(walk.getElement());
            walk = walk.getNext();
            if (walk != null)
                sb.append(", ");
        }
        sb.append(")");
        return sb.toString();
    }

}

// Given in programming assignment 1.
class Node<E> {
    E element; // element
    Node<E> prev; // preceeding Node
    Node<E> next; // succeeding Node
    Cluster<E> cluster;

    Node(Node<E> prev, E element, Node<E> next, Cluster cluster) {
        this.prev = prev;
        this.element = element;
        this.next = next;
        this.cluster = cluster;
    }

    // getters
    public E getElement() {
        return this.element;
    }

    public Cluster<E> getCluster() {
        return this.cluster;
    }

    public Node<E> getNext() {
        return next;
    }

    public Node<E> getPrev() {
        return prev;
    }

    // setters
    public void setElement(E elem) {
        this.element = elem;
    }

    public void setCluster(Cluster<E> n) {
        this.cluster = n;
    }

    public void setNext(Node<E> n) {
        next = n;
    }

    public void setPrev(Node<E> n) {
        prev = n;
    }

    @Override
    public String toString() {
        return "Node(" + this.element + ")";
    }

}

// Given in programming assignment 1.
class Partition<E> {
    /**
     * Implements Partition as a linked list
     */

    Cluster<E> pheader;
    Cluster<E> ptrailer;
    int size = 0;

    public Partition() {
        pheader = new Cluster<>(null, null);
        ptrailer = new Cluster<>(pheader, null);
        pheader.setNextCluster(ptrailer);

    }

    public Node<E> makeCluster(E x) {
        Cluster<E> c = new Cluster<>(ptrailer.getPrevCluster(), ptrailer);
        Node<E> pos = c.addLast(x);
        this.addLast(c);
        return pos;

    }

    public Node<E> find(Node<E> n) {

        return n.getCluster().first();

    }

    public void union(Node<E> p, Node<E> q) {

        Cluster<E> a = find(p).getCluster();
        Cluster<E> b = find(q).getCluster();
        if (a == b) {
            return;
        }
        Cluster<E> c1 = min(a, b); // cluster to get merged
        Cluster<E> c2 = c1.equals(a) ? b : a; // cluster to merge into

        c2.add(c1);

        Node<E> curr = c2.first();
        while (curr != null) {
            curr.setCluster(c2);
            curr = c2.after(curr);
        }

        remove(c1);

    }

    public E element(Node<E> n) {
        return n.getElement();

    }

    public int clusterSize(Node<E> n) {
        Node<E> leader = find(n);
        return leader.getCluster().size();

    }

    public Node<E>[] clusterPositions(Node<E> p) {
        Node<E> n = find(p);
        Cluster<E> c = n.getCluster();
        Node<E>[] positions = new Node[c.size()];
        Node<E> curr = c.first();
        for (int i = 0; i < c.size(); i++) {
            positions[i] = curr;
            curr = c.after(curr);

        }
        return positions;

    }

    public Integer[] clusterSizes() {
        Cluster<E> curr = pheader.getNextCluster();
        Integer[] sizes = new Integer[size()];
        for (int i = 0; i < size(); i++) {
            sizes[i] = curr.size();
            curr = after(curr);

        }
        Arrays.sort(sizes, Collections.reverseOrder());

        return sizes;

    }

    public int numberOfClusters() {
        return size;
    }

    // Utility methods
    public int size() {
        return size;
    }

    public void addLast(Cluster<E> c) {
        ptrailer.getPrevCluster().setNextCluster(c);
        ptrailer.setPrevCluster(c);
        size++;

    }

    public Cluster<E> after(Cluster<E> c) {
        return c.getNextCluster();

    }

    public Cluster<E> before(Cluster<E> c) {
        return c.getPrevCluster();

    }

    public Cluster<E> first() {
        return pheader.getNextCluster();
    }

    public void remove(Cluster<E> c) {
        c.getPrevCluster().setNextCluster(c.getNextCluster());
        c.getNextCluster().setPrevCluster(c.getPrevCluster());
        size--;
    }

    // Helper methods
    public Cluster<E> min(Cluster<E> c1, Cluster<E> c2) {
        return (c1.size() <= c2.size() ? c1 : c2);
    }

    public Node<E>[] getLeaders() {
        Node<E>[] leaders = new Node[size()];
        Cluster<E> curr = first();
        for (int i = 0; i < size(); i++) {
            leaders[i] = curr.first();
            curr = after(curr);
        }
        return leaders;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("Partition\n|\n");
        Cluster<E> curr = this.first();
        while (curr != ptrailer) {
            sb.append(curr.toString());
            curr = this.after(curr);
            if (curr != ptrailer)
                sb.append("\n");
        }
        return sb.toString();
    }

}

// Implementing Edge to keep track of edges. Implementing Comparable to allow
// for edges to be sorted in order.
class Edge<E> implements Comparable<Edge<E>> {

    // Required attributes
    private Node<E> u;
    private Node<E> v;
    private int weight;

    // Constructor
    public Edge(Node<E> u, Node<E> v, int weight) {
        this.u = u;
        this.v = v;
        this.weight = weight;
    }

    // Getting verticie 1
    public Node<E> getU() {
        return u;
    }

    // Getting verticie 2
    public Node<E> getV() {
        return v;
    }

    // Getting weight
    public int getWeight() {
        return weight;
    }

    // Implementing compareto
    @Override
    public int compareTo(Edge<E> other) {
        return Integer.compare(this.weight, other.weight);
    }

    // Implementing toString (used in assignment as debugger)
    @Override
    public String toString() {
        return "(" + u.getElement() + " - " + v.getElement() + ", " + weight + ")";
    }
}

/*
 * Helper classes for the metro graph and Dijkstra algorithm.
 */

class MetroEdge {
    int to;
    int weight;

    MetroEdge(int to, int weight) {
        this.to = to;
        this.weight = weight;
    }
}

class NodeDist {
    int v;
    int dist;

    NodeDist(int v, int dist) {
        this.v = v;
        this.dist = dist;
    }
}
