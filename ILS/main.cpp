#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>

/*
 typedef struct Solucao{
    vector<int> sequence;
    double valorObj;
 } Solucao;
Assim que é a struct de solução em C*/

/*
    CL = V\V' significa que CL é igual ao conjunto de vértices de V que não estão em V'
    CL, então, seria o complemento de V'
*/

/*
    Falta ajeitar o update solution
*/

const int numberOfVertices = 10;
std::vector<int> vertices;
int edgeCost[numberOfVertices][numberOfVertices];

class Solution{
    private:
        std::vector<int> sequence;
        double cost;

    public:

        Solution(std::vector<int> seq, double obj):sequence{seq}, cost{obj}{}

        Solution(){}

        void print(){
            for(int i = 0; i < sequence.size(); i++){
                std::cout << sequence.at(i);
                if(i != sequence.size() - 1){
                    std::cout << " -> ";
                }
            }
            std::cout << std::endl;
        }

        void setSequence(std::vector<int> seq){
            sequence = seq;
        }

        std::vector<int> getSequence(){
            return sequence;
        }

        void setCost(double obj){
            cost = obj;
        }

        double getCost(){
            return cost;
        }

        void pushBackSequence(int value){
            sequence.push_back(value);
        }
};

class InsertionInfo{
    private:
        int insertedNode;
        int removedEdge;
        double cost;
    
    public:
        int getInsertedNode(){
            return insertedNode;
        }

        int getRemovedEdge(){
            return removedEdge;
        }

        double getCost(){
            return cost;
        }

        void setInsertedNode(int node){
            insertedNode = node;
        }

        void setRemovedEdge(int edge){
            removedEdge = edge;
        }

        void setCost(double value){
            cost = value;
        }

        void print(){
            printf("{ insertedNode: %d, removedEdge: %d, cost: %f }\n", insertedNode, removedEdge, cost);
        }
};

Solution iteratedLocalSearch(int maxIter, int maxIterILS);

Solution construction();

Solution perturbation(Solution bestSolution);

void localSearch(Solution* possibleSolution);

bool greaterCost(Solution first, Solution second);

void whileILS(Solution solution, int *iterILS);

std::vector<InsertionInfo> insertionCostCalculation(Solution &solution, std::vector<int> &complement);

std::vector<int> choose3RandomNodes();

std::vector<int> remainingNodes(std::vector<int> sequence);

void insertionSort(std::vector<InsertionInfo> &sequence);

void updateSolution(Solution &solution, int choosen, std::vector<int> &complement);


int main(void){

    for(int i = 0; i < numberOfVertices; i++){
        vertices.push_back(i + 1);
    }

    Solution solution = construction();

    printf("\n");
    solution.print();

    // Solution solution;
    // solution.pushBackSequence(vertices.at(0));
    // solution.pushBackSequence(vertices.at(1));
    // solution.pushBackSequence(vertices.at(2));

    // std::vector<int> complement;
    // for(int i = solution.getSequence().size(); i < numberOfVertices; i++){
    //     complement.push_back(i);
    // }

    // std::vector<InsertionInfo> insertionCost = insertionCostCalculation(solution, complement);

    // for(auto insertion: insertionCost){
    //     insertion.print();
    // }

    return 0;
}

Solution construction(){
    Solution solution;

    solution.setSequence(choose3RandomNodes());
    std::vector<int> complement = remainingNodes(solution.getSequence());

    double alpha;
    int choosen;

    while(!complement.empty()){
        std::vector<InsertionInfo> insertionCost = insertionCostCalculation(solution, complement);
        insertionSort(insertionCost);
        alpha = (double) rand() / RAND_MAX;
        choosen = rand() % ((int) ceil (alpha * insertionCost.size()));
        updateSolution(solution, insertionCost.at(choosen), complement);
    }

    return solution;
}

Solution perturbation(Solution bestSolution){
    Solution solution;
    return solution;
}

void localSearch(Solution* possibleSolution){

}

bool greaterCost(Solution first, Solution second){
    return first.getCost() < second.getCost();
}

Solution iteratedLocalSearch(int maxIter, int maxIterILS){
    Solution bestOfAll;
    int iterILS;

    bestOfAll.setCost(INFINITY);

    for(int i = 0; i < maxIter; i++){
        Solution solution = construction(); //Construção de uma solução baseado em "palpites educativos"
        Solution best = solution;

        iterILS = 0;

        while(iterILS <= maxIterILS){ //Em tese, a local search em si fica aqui
            localSearch(&solution); //Tentando fazer melhoras na solução através de pequenas modificações

            if(greaterCost(solution, best)){
                best = solution;
                iterILS = 0; //Note que o iterILS vira zero quando uma solução melhor é encontrada
            }

            solution = perturbation(best); //Provoca uma "pertubação" na solução para melhorá-la
            iterILS++; 
        }

        if(greaterCost(best, bestOfAll)){
            bestOfAll = best;
        }
    }

    return bestOfAll;
}

std::vector<InsertionInfo> insertionCostCalculation(Solution &solution, std::vector<int> &complement){
    int solutionSequenceSize = solution.getSequence().size();

    std::vector<InsertionInfo> insertionCost( (solutionSequenceSize) * (complement.size()));

    int l = 0;
    for(int a = 0, b = 1; a < solutionSequenceSize - 1; a++, b++){
        int vertexI = solution.getSequence().at(a);
        int vertexJ = solution.getSequence().at(b);

        for(auto vertexK : complement){
            insertionCost[l].setCost(edgeCost[vertexI][vertexK] + edgeCost[vertexJ][vertexK] - edgeCost[vertexI][vertexJ]);
            insertionCost[l].setInsertedNode(vertexK);
            insertionCost[l].setRemovedEdge(a);
            l++;
        }
    }

    return insertionCost;
}

std::vector<int> choose3RandomNodes(){
    std::vector<int> ret;

    for(int i = 0; i < 3; i++){
        ret.push_back( rand() % numberOfVertices);
    }

    return ret;
}

std::vector<int> remainingNodes(std::vector<int> sequence){
    std::vector<int> ret;

    for(int i = 0; i < numberOfVertices; i++){
        if(std::find(sequence.begin(), sequence.end(), i) == sequence.end()){
            ret.push_back(i);
        }
    }

    return ret;
}

void insertionSort(std::vector<InsertionInfo> &sequence){
    int j;
    InsertionInfo temp;
    for(int i = 1; i < sequence.size(); i++){
        temp = sequence.at(i);
        j = i - 1;
        while(j >= 0 && sequence.at(j).getCost() < temp.getCost()){
            sequence[j + 1] = sequence[j];
            j--;
        }
        sequence[j + 1] = temp;
    }
}

void updateSolution(Solution &solution, InsertionInfo choosen, std::vector<int> &complement){
    solution.pushBackSequence(choosen);

    for(int i = 0; i < complement.size(); i++){
        if(complement.at(i) == choosen){
            complement.erase(complement.begin() + i);
            break;
        }
    }

}





