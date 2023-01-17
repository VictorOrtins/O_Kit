#include <iostream>
#include <vector>

class Solution{
    private:
        std::vector<int> sequencia;
        double valorObj;

    public:

        Solution(std::vector<int> seq, double obj):sequencia{seq}, valorObj{obj}{}

        void exibirSolucao(){
            for(int i = 0; i < sequencia.size(); i++){
                std::cout << sequencia.at(i);
                if(i != sequencia.size() - 1){
                    std::cout << " -> ";
                }
            }
            std::cout << std::endl;
        }

        void setSequencia(std::vector<int> seq){
            sequencia = seq;
        }

        std::vector<int> getSequencia(){
            return sequencia;
        }

        void setValorObj(double obj){
            valorObj = obj;
        }

        double getValorObj(){
            return valorObj;
        }

        void pushBackSequencia(int value){
            sequencia.push_back(value);
        }
}; 

int main(void){


    return 0;
}