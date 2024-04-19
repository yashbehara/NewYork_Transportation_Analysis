/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package newyork_transportation;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;

import java.util.*;

public class NewYork_Transportation {
    
    private static int currentDatasetChoice;
    public static void main(String[] args) {
        int datasetChoice = getUserDatasetChoice();
        currentDatasetChoice = datasetChoice;
        Graph graph = loadGraphFromDataset(datasetChoice);
        Scanner scanner = new Scanner(System.in);

        while (true) {
            int planChoice = getUserPlanChoice();
            switch (planChoice) {
                case 1:
                    handleSingleStation(graph);
                    break;
                case 2:
                    handleOptimalRoute(graph);
                    break;
                case 3:
                    handleCompleteNetwork(graph);
                    break;
                default:
                    System.out.println("Invalid plan choice. Please try again.");
            }
            System.out.println("\n");
            System.out.println("Do you want to continue? (yes/no)");
            String continueChoice = scanner.next();
            if (!continueChoice.equalsIgnoreCase("yes")) {
                break;
            }
        }
    }
    
    public static Graph loadGraphFromDataset(int datasetChoice) {
        String csvFile = "";
        switch (datasetChoice) {
            case 1:
                csvFile = "/Users/yash/NetBeansProjects/NewYork_Transportation/src/newyork_transportation/Cost_Dataset.csv";
                break;
            case 2:
                csvFile = "/Users/yash/NetBeansProjects/NewYork_Transportation/src/newyork_transportation/Distance_Dataset.csv";
                break;
            case 3:
                csvFile = "/Users/yash/NetBeansProjects/NewYork_Transportation/src/newyork_transportation/Time_Dataset.csv";
                break;
            case 4:
                csvFile = "/Users/yash/NetBeansProjects/NewYork_Transportation/src/newyork_transportation/CTD_Dataset.csv";
                break;
            case 5:
                csvFile = "/Users/yash/NetBeansProjects/NewYork_Transportation/src/newyork_transportation/DTC_Dataset.csv";
                break;
            case 6:
                csvFile = "/Users/yash/NetBeansProjects/NewYork_Transportation/src/newyork_transportation/TCD_Dataset.csv";
                break;
            default:
                throw new IllegalArgumentException("Invalid dataset choice.");
        }

        Graph graph = new Graph();
        try (BufferedReader br = new BufferedReader(new FileReader(csvFile))) {
            String line;
            boolean firstLine = true;
            String[] headers = null;

            while ((line = br.readLine()) != null) {
                String[] values = line.split(",");
                if (firstLine) {
                    headers = new String[values.length - 1];
                    System.arraycopy(values, 1, headers, 0, values.length - 1);
                    firstLine = false;
                } else {
                    String street = values[0];
                    for (int i = 1; i < values.length; i++) {
                        String station = headers[i - 1];
                        double weight = Double.parseDouble(values[i]);
                        graph.addEdge(street, station, weight);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return graph;
    }
    
    public static int getUserPlanChoice() {
            Scanner scanner = new Scanner(System.in);
            System.out.println("\nChoose a plan (Single Station : 1, Optimal Route : 2, Complete Network : 3):");
            return scanner.nextInt();
    }
    
    
    public static void handleSingleStation(Graph graph) {
        Scanner scanner = new Scanner(System.in);

        // List all stations
        System.out.println("List of stations:");
        Set<String> stationsList = graph.getVertices();
        for (String station : stationsList) {
            System.out.println(station);
        }

        // Prompt user to choose source station
        String sourceStreet;
        boolean validStation = false;
        do {
            System.out.println("\nChoose a source station:");
            sourceStreet = scanner.nextLine();
            if (stationsList.contains(sourceStreet)) {
                validStation = true;
            } else {
                System.out.println("Please enter a valid station.");
            }
        } while (!validStation);

        // Prompt user to choose analysis type
        int analysisChoice;
        do {
            System.out.println("\nChoose analysis type:");
            System.out.println("1. Direct Analysis");
            System.out.println("2. Balanced Analysis of Cost, Distance, and Time");
            analysisChoice = scanner.nextInt();
            if (analysisChoice != 1 && analysisChoice != 2) {
                System.out.println("Invalid choice. Please enter 1 or 2.");
            }
        } while (analysisChoice != 1 && analysisChoice != 2);

        if (analysisChoice == 1) {
            // Algorithm 1: Dijkstra's Algorithm
            Map<String, Double> shortestDistances = dijkstra(graph, sourceStreet);
            System.out.println("\nShortest distances from " + sourceStreet + ":");
            for (Map.Entry<String, Double> entry : shortestDistances.entrySet()) {
                System.out.println(entry.getKey() + ": " + entry.getValue());
            }
        } else if (analysisChoice == 2) {
            // Algorithm 2: Balanced Analysis of Cost, Distance, and Time
            Map<String, Double> averageShortestDistances = performBalancedAnalysis(graph, sourceStreet);
            System.out.println("\nAverage shortest distances from " + sourceStreet + ":");
            for (Map.Entry<String, Double> entry : averageShortestDistances.entrySet()) {
                System.out.println(entry.getKey() + ": " + entry.getValue());
            }
        }
    }
    //    Method to perform Balanced analysis ( avg of cost, time and distance ) on Dijkstra's algorithm
    public static Map<String, Double> performBalancedAnalysis(Graph graph, String sourceStreet) {
        Map<String, Double> averageShortestDistances = new HashMap<>();

        // Perform Dijkstra's algorithm for each dataset and calculate average shortest distances
        for (int i = 1; i <= 3; i++) {
            Graph datasetGraph = loadGraphFromDataset(i);
            Map<String, Double> shortestDistances = dijkstra(datasetGraph, sourceStreet);

            for (Map.Entry<String, Double> entry : shortestDistances.entrySet()) {
                String station = entry.getKey();
                double distance = entry.getValue();
                averageShortestDistances.merge(station, distance, Double::sum);
            }
        }

        // Calculate average distances and round to two decimal points
        DecimalFormat df = new DecimalFormat("#.##");
        for (Map.Entry<String, Double> entry : averageShortestDistances.entrySet()) {
            double averageDistance = entry.getValue() / 3; // Dividing by 3 since we analyzed 3 datasets
            String roundedDistance = df.format(averageDistance);
            entry.setValue(Double.parseDouble(roundedDistance));
        }

        return averageShortestDistances;
    }

    public static void handleOptimalRoute(Graph graph) {
            Set<Edge> mstEdges = kruskalMST(graph);
            System.out.println("Minimum cost to cover all the stations:");
            System.out.println("\n");

            for (Edge edge : mstEdges) {
                System.out.println(edge);
            }
            
            // Perform analysis
            performAnalysis(mstEdges);
        }

    public static void handleCompleteNetwork(Graph graph) {
            int[][] shortestDistancesMatrix = floydWarshall(graph);
            System.out.println("Shortest distances between every pair of stations:");
            Set<String> stations = graph.getVertices();
            for (int i = 0; i < shortestDistancesMatrix.length; i++) {
                System.out.println("\nFrom " + stations.toArray()[i]);
                for (int j = 0; j < shortestDistancesMatrix[i].length; j++) {
                    if (i != j) {
                        System.out.println("To " + stations.toArray()[j] + ": " + shortestDistancesMatrix[i][j]);
                    }
                }
            }
    }


    public static int getUserDatasetChoice() {
            Scanner scanner = new Scanner(System.in);
            System.out.println("Choose a dataset (Cost : 1, Distance : 2, Time : 3, CTD (Cost, Time, Distance): 4, DTC (Distance, Time, Cost): 5, TCD (Time, Cost, Distance): 6");
            return scanner.nextInt();
    }

    public static Map<String, Double> dijkstra(Graph graph, String source) {
        Map<String, Double> distances = new HashMap<>();
        Set<String> visited = new HashSet<>();
        PriorityQueue<Node> pq = new PriorityQueue<>(Comparator.comparingDouble(node -> node.distance));
        distances.put(source, 0.0);
        pq.offer(new Node(source, 0.0));

        while (!pq.isEmpty()) {
            Node current = pq.poll();
            String currentVertex = current.vertex;
            if (visited.contains(currentVertex)) {
                continue;
            }
            visited.add(currentVertex);

            Map<String, Double> neighbors = graph.adjacencyList.get(currentVertex);
            if (neighbors != null) {
                for (Map.Entry<String, Double> neighborEntry : neighbors.entrySet()) {
                    String neighbor = neighborEntry.getKey();
                    double weight = neighborEntry.getValue();

                    // Skip neighbors with weight 0
                    if (weight == 0.0) {
                        continue;
                    }
                    double newDistance = distances.get(currentVertex) + weight;
                    if (!distances.containsKey(neighbor) || newDistance < distances.get(neighbor)) {
                        distances.put(neighbor, newDistance);
                        pq.offer(new Node(neighbor, newDistance));
                    }
                }
            }
        }

        return distances;
    }

    public static Set<Edge> kruskalMST(Graph graph) {
        Set<Edge> mst = new HashSet<>();
        PriorityQueue<Edge> pq = new PriorityQueue<>(Comparator.comparingDouble(edge -> edge.weight));

        // Add all edges to priority queue
        for (String vertex : graph.getVertices()) {
            Map<String, Double> neighbors = graph.adjacencyList.get(vertex);
            if (neighbors != null) {
                for (Map.Entry<String, Double> neighborEntry : neighbors.entrySet()) {
                    String neighbor = neighborEntry.getKey();
                    double weight = neighborEntry.getValue();
                    if (weight != 0.0) {
                        pq.offer(new Edge(vertex, neighbor, weight));
                    }
                }
            }
        }

        DisjointSet disjointSet = new DisjointSet();
        for (String vertex : graph.getVertices()) {
            disjointSet.makeSet(vertex);
        }

        while (!pq.isEmpty() && mst.size() < graph.getVertices().size() - 1) {
            Edge edge = pq.poll();
            String source = edge.source;
            String destination = edge.destination;
            if (disjointSet.findSet(source) != disjointSet.findSet(destination)) {
                mst.add(edge);
                disjointSet.union(source, destination);
            }
        }

        return mst;
    }

    public static int[][] floydWarshall(Graph graph) {
        int n = graph.getVertices().size();
        int[][] distances = new int[n][n];

        // Initialize distances array with infinity
        int infinity = Integer.MAX_VALUE / 2; // Avoid overflow
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                distances[i][j] = infinity;
            }
        }

        // Fill in initial distances based on the graph
        for (String source : graph.getVertices()) {
            int sourceIndex = getIndex(graph.getVertices(), source);
            Map<String, Double> neighbors = graph.adjacencyList.get(source);
            if (neighbors != null) {
                for (Map.Entry<String, Double> neighborEntry : neighbors.entrySet()) {
                    String destination = neighborEntry.getKey();
                    double weight = neighborEntry.getValue();
                    if (weight != 0.0) { // Skip 0 weights
                        int destinationIndex = getIndex(graph.getVertices(), destination);
                        distances[sourceIndex][destinationIndex] = (int) weight;
                    }
                }
            }
        }

        // Apply Floyd-Warshall algorithm
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (distances[i][k] != infinity && distances[k][j] != infinity &&
                            distances[i][k] + distances[k][j] < distances[i][j]) {
                        distances[i][j] = distances[i][k] + distances[k][j];
                    }
                }
            }
        }

        return distances;
    }
    
    public static int getIndex(Set<String> vertices, String vertex) {
        int index = 0;
        for (String v : vertices) {
            if (v.equals(vertex)) {
                return index;
            }
            index++;
        }
        return -1;
    }

    static class Graph {
        private Map<String, Map<String, Double>> adjacencyList;

        public Graph() {
            adjacencyList = new HashMap<>();
        }

        public void addEdge(String source, String destination, double weight) {
            adjacencyList.putIfAbsent(source, new HashMap<>());
            adjacencyList.get(source).put(destination, weight);
            adjacencyList.putIfAbsent(destination, new HashMap<>());
        }

        public Set<String> getVertices() {
            return adjacencyList.keySet();
        }
    }

    static class Node {
        String vertex;
        double distance;

        public Node(String vertex, double distance) {
            this.vertex = vertex;
            this.distance = distance;
        }
    }

    static class Edge {
        String source;
        String destination;
        double weight;

        public Edge(String source, String destination, double weight) {
            this.source = source;
            this.destination = destination;
            this.weight = weight;
        }

        @Override
        public String toString() {
            return "(" + source + " -> " + destination + ": " + weight + ")";
        }
    }

    static class DisjointSet {
        private Map<String, String> parent;

        public DisjointSet() {
            parent = new HashMap<>();
        }

        public void makeSet(String vertex) {
            parent.put(vertex, vertex);
        }

        public String findSet(String vertex) {
            if (!parent.get(vertex).equals(vertex)) {
                parent.put(vertex, findSet(parent.get(vertex)));
            }
            return parent.get(vertex);
        }

        public void union(String vertex1, String vertex2) {
            String root1 = findSet(vertex1);
            String root2 = findSet(vertex2);
            if (!root1.equals(root2)) {
                parent.put(root1, root2);
            }
        }
    }
    
    public static void performAnalysis(Set<Edge> mstEdges) {
        // Total cost, distance, or time
        double totalValue = 0.0;

        // Longest and shortest edge
        Edge longestEdge = null;
        Edge shortestEdge = null;

        // Station connectivity
        Map<String, Integer> stationConnectivity = new HashMap<>();

        // Initialize average edge weight
        double totalEdgeWeight = 0.0;
        int edgeCount = 0;

        // Iterate through the edges in the MST
        for (Edge edge : mstEdges) {
            // Update total value (cost, distance, or time)
            totalValue += edge.weight;

            // Update longest and shortest edge
            if (longestEdge == null || edge.weight > longestEdge.weight) {
                longestEdge = edge;
            }
            if (shortestEdge == null || edge.weight < shortestEdge.weight) {
                shortestEdge = edge;
            }

            // Update station connectivity
            stationConnectivity.put(edge.source, stationConnectivity.getOrDefault(edge.source, 0) + 1);
            stationConnectivity.put(edge.destination, stationConnectivity.getOrDefault(edge.destination, 0) + 1);

            // Update average edge weight
            totalEdgeWeight += edge.weight;
            edgeCount++;
        }

        // Calculate average edge weight
        double averageEdgeWeight = totalEdgeWeight / edgeCount;

        // Display analysis results
        System.out.println("\nAnalysis Results:");
        
        System.out.println("\n1. Total Cost of the path: " + totalValue);


        System.out.println("\n2. Longest Edge: " + longestEdge);
        System.out.println("   Shortest Edge: " + shortestEdge);
        System.out.println("\n3. Average Edge Weight: " + averageEdgeWeight);
        System.out.println("\n4. Station Connectivity:");
        for (Map.Entry<String, Integer> entry : stationConnectivity.entrySet()) {
            System.out.println("   " + entry.getKey() + ": " + entry.getValue() + " connections");
        }
    }

}