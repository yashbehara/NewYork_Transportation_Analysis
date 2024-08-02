# NewYork_Transportation_Analysis


## Summary : 

1. The objective of this project is to develop a transportation network that analyses routes between stations, considering distance, cost, and time constraints. By constructing a transportation network graph and applying graph algorithms, the project aims to identify optimal routes that minimize transportation costs while meeting specified time and distance requirements.

2. This project gives user the flexibility to choose the weights of the graph. It can be Distance, Cost or Time. Taking user’s input into consideration, the dataset will be loaded and the graph will be developed.

## Implementation : 
1.	Dataset is about NewYork. 
    a.	Vertexes are Stations
    b.	Edges are routes between them
    c.	Weights:
        i.	Cost
        ii.	Distance
        iii.	Time
2.	There are totally 6 different datasets and user will be given preference to choose.
    a.	Cost based weights
    b.	Distance based
    c.	Time based
    d.	CDT (Cost Distance Time)
    e.	DTC (Distance Time Cost) 
    f.	TCD (Time Cost Distance)
3.	Based on the constraint chosen, the analysis will be performed


## Project Work Flow : 

1. User after selecting the dataset ( 1. Cost, 2. Distance, 3. Time ): 
	a. Selects Single Station ( For shortest path )
		i. After entering Source Station
			1. If Direct Analysis is chosen
				a. Direct Dijkstra’s algorithm is performed for shortest paths 
			2. If Balanced Analysis is chosen
				a. Average of 3 datasets (cost, distance, time) is performed and then shortest paths are displayed
	b. Selects Optimal Path ( Minimum Spanning Tree )
		i. Analysis is performed after Kruskal’s algorithm is executed
			1. Total cost of path
			2. Longest & Shortest edges
			3. Avg weight of the edge
			4. No. of connections of each station
	c. Selects Complete Network ( All sources shortest paths )
		i. Analysis is performed with Floyd Warshall’s algorithm
			1. Distance from one stations to all other stations is shown
			2. This is repeated for all possible sources in the network
	d. Selects DFS ( Allows user to enter stations they wish to avoid )
		i. Enter Starting Station
		ii. Enter Station to avoid
		iii. Analysis is performed with all stations traversed except the avoiding station
	e. Selects BFS
		i. Full network traversal is performed


## Algorithms Implemented : 

1. Dijkstra’s Algorithm
2. Kruskal’ s Algorithm
3. Floyd Warshall Algorithm
4. Depth First Search
5. Breadth First Search
