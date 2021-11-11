# Steady-State-Network-Flow
This gives steady state Network Flow from the given Pref-Defined Network Flow of a System.


## Results of Random Testing
* 100 Random networks with 100 Random Pages weight
```python
Comparing with PageRank Algorithm
max percentage error [0.15969399]  min percentage error  [-0.1725098]
Comparing with Power Iteration Algorithm
max percentage error  [0.15969399]  min percentage error  [-0.1725098]

±0.18% error found

```
### What is Network Flow?

* A network is a directed graph that forms a system of nodes and edges that connect them to each other. 
* A network has multiple connotations across Graph Theory, Linear Algebra, Electric flow, Communication. 
* Network flow is essentially a graph theory problem where it defines the capacity of each edge in a connected network to transit data

### What is Steady State in Network Flow?

* The state of a system in which the flow of network does not affect the population of nodes.
* The stabilized flow of the network results in such a flow that the relative population at the nodes are plateaued according to the potential. 
* Which implies that the number of visitors / data at the nodes are unchanged irrespective of the flow after a certain number of iterations.


### What the result describes?

* We have defined our primary subject matter as steady state analysis of a network flow. 
* Steady state analysis means repetitive analysis of a flow in network till it reaches the steady state.
* As in PageRank we entertain random clicks/visits to a particular link/node of network or with the real time data sets we can get "more important page" as more scored PageRank and less important page as less scored PageRank. 
* Whereas, in network flow we have some predefined flow, after over and over iterations we come across a steady state which is steady state of a network flow. 
* In steady state the value of network flow to the node signifies its influence in the network which can be associated to the PageRank as a rank in the network.

### Assumptions and Limitations: -

* The main assumption we use is that since we use Jacobi Algorithm to find the eigenvectors and eigen-values we work with Symmetric Matrices only (By using Python Libraries we can eradicate this limitation).
* The test cases are preferably in range of 2 x 2 to 5 x 5 matrices as we use matrix multiplication algorithm in O(n^3) time complexity which uses 3 For loops(We can make it more efficient and optimized by using libabries).
* The sum of elements of each column of input matrix must be Unity. Each element input must be non negative (It is the assumption we can sacle our Code accordingly).

### Approach for calculation of Steady state : -

* We have used python to code the Network flow and it’s steady state flow.
* We find out the eigenvalues  and eigen vectors of the matrix
* Now we decompose the matrix to  (S)*(D)*(𝑺^(−𝟏)) where  S contains the eigen vectors of A; Whereas D is diagonal matrix containing eigenvalues of A
* Now we scale it to kth degree so that 〖(𝑨〗^𝒌)=(𝑺𝑫^𝒌 𝑺^(−𝟏))
* Now once we get 〖(𝑨〗^𝒌⋅𝑼_𝟎)=〖(𝑼〗_𝒌) which was our equation 
* After a certain number of iterations when the value starts to plateau we can say that the network has attained a steady state.
* This means the vector that we finally get is the vector that represents the rank of the nodes in the given graph, which is useful in ranking them from the highest to lowest.

### How the code works? 

* We input the input matrix in the code with respect to the assumptions and limitations. 
* Then we find the eigenvalues and eigenvectors of the input matrix. 
* The eigen space then is further diagonalized with ‘k’ iterations up until it reaches steady state.
* Now here in our strategy to reduce time complexity we apply the fact that while scaling (S)*(D)*(𝑺^(−𝟏)) and multiplying it again with (S)*(D)*(𝑺^(−𝟏)) the (S)(𝑺^(−𝟏)) results  for an identity matrix. 
* This optimizes our approach as we use the naïve multiplication algorithm which takes O(n^3) complexity

### Linear Algebra Concepts used for Project: -

* In order to find Eigen values and Eigen vectors  we use Jacobi algorithm.
* Then we decompose the matrix  in terms of its diagonal matrix and eigenvectors.
* For the above , we mainly define two functions diagonalization  to find the  product of eigenvectors matrix, Diagonal matrix and inverse of eigenvectors matrix) and Inverse function to find inverse of eigen vectors matrix.
* In Inverse function, we find the inverse of matrix using elementary  Row Operations.

### Why this approach works?

* This approach works because the probability of number of visitors/ data going from one node to another depends upon the previous iteration. 
* The flow also depends upon a node's connection with other important nodes 
* This approach works because of the. In a matrix, the probability is are Eigen Centrality Concept designated to nodes. There will be a dominant Eigen Value for which the flow attains a Steady State
* We capitalize on this fact that rank of a node only depends on the previous rank. So we find the eigen values and eigenvectors to exponentiation in order to iterate.

### Inferences from the Code: -

* The symmetric matrices do give a steady state where the vectors that represents the Steady State is (1, 1, 1, ….). 
* The non- symmetric matrices do give out the percentage of population/ data that is present at the saturation/ steady state of the flow. 
* The steady state also speaks about page rank in which the value of each node which is a webpage shows what numbers of visitors are going to stay on a given node and what number are going to transit. 
* This helps in deciding the relevancy of the pages.

## Conclusion: -

* The Steady State of a network is indeed an (proportional) eigen vector of the input Matrix
* The Page Rank Algorithm and Power Iteration Algorithm give out similar results as our code because of Eigen Centrality.
* Repetitive analysis does actually give out the steady state of the network flow and so does it give out convincing results for algorithms like PageRank.

