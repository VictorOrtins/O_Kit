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

    Solution s1{{1,6,3,2,5,4,1}, 0};
    s1.exibirSolucao();

    Solution s2{{}, 0.0};
    s2.pushBackSequencia(1);
    s2.pushBackSequencia(3);
    s2.pushBackSequencia(1);

    s2.exibirSolucao();

    //Tá emperrando nisso aqui
    s2.getSequencia().insert(s2.getSequencia().begin() + 1, 4);
    s2.getSequencia().insert(s2.getSequencia().end() - 1, 5);
    s2.getSequencia().insert(s2.getSequencia().begin() + 2, 2);

    s2.exibirSolucao();


    return 0;
}