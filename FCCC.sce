//----------------------------------------------
// Algorítimo OCF responsável pela otimização na
// corrente injetada na rede distribuição.
// Projeto de Pesquisa FAPESP
// Projeto número: #2019/24128-2
// @date 01/07/2020
// @author Higor de Paula Kolecha
// @author Adolfo Blengini Neto
// @author Marcius Fabius Henriques de Carvalho.
// @version 1.0
//----------------------------------------------

// Responsável pela limpeza toda memória
clear;
// Responsável pela limpeza de tela
clc;

//Solicitação ao usuário do endereço para obtenção do arquivo de entrada.
entradaDeDados=input("Digite o endereço com localização do arquivo de entrada de dados. Obs: seguir instruções no arquivo Guideline. ");

// Arquivo de Entrada com a estrutura da Rede de Distribuição.
M = fscanfMat(entradaDeDados, "%lg"); // Importa arquivo que apresenta os dados de entrada da rede. 

//"C:\Users\Higor\Desktop\rede5barras.txt"
//"C:\Users\Higor\Desktop\Iniciação científica\IC 2020 marcius scilab\Programação\11 barras v2\rede11barras.txt"
//"C:\Users\Higor\Documents\OtimizacaoDeRedes\Arquivos de Entrada\FCCC\rede34barras2.txt"
// Salva o número de ramos
NR=M(1,1); 
// Salva o número de barras
NB=M(1,2); 
//valor de refencia de tensão inicial
vf=1.0;
//valor de convergencia
e=0.0001;
// -----------------------------------------------------
// Funcao responsável pela construçao do laço externo 
// Entrada : Estrutura do arquivo completo 
// Saida : Matriz de lacos externos
// -----------------------------------------------------

function [quaislacos, barralaco]=lacoE(M)
    //Entrada
    M(1,1)=1; // Inicia os contadores de linhas e colunas
    M(1,2)=1; // Inicia os contadores de linhas e colunas
    
    O=M(:,1);  // Criação do vetor origem.
    D=M(:,2);  // Criação do vetor destino.
    
    quaislacos=[] // Criação do caminho que a corrente faz para poder criar um laço.
    
    barralaco=[];
    for i=NR:-1:1
        buss=i;     
        if (M(i,3)~=0) //Verificação da existencia de carga. Caso afirmativo, passa para a criação dos laços
            laco=zeros(NR,1); // Criação do caminho que a corrente faz para poder criar um laço.
            barralaco=cat(1,barralaco,i);// Faz a concatenação da matriz "barralaco" a fim de serem adicionadas apenas laços existentes

            //Busca amontante até feeder, ou seja, barra 1
            while (buss~=1)
                o=find(D(:,1)==buss); //Mostra a posição onde se encontra o próximo ramo para formação de laços
                laco(o(1,1),1)=M(o(1,1),5); //Monta o vetor caminhos, mostrando os laços feitos
                buss=O(o(1,1)); //Recebe o valor que se encontra na posição acima destacada dentro do vetor origem
            end
            quaislacos = cat(2,quaislacos,laco);
        end
    end    
    l = size(quaislacos,"c")// Mede o tamanho da matriz quaislaacos
    complemento = eye(l,l);
    quaislacos = cat(1,quaislacos, complemento)// Faz a concatenação da matriz quais lacos com a matriz complemento
endfunction //Fim da função de criação de laços

// ---------------------------------------------
// Funcao responsavel pela criação da Matriz A
//----------------------------------------------

function MA=MatrizA(M, LE)
    n = size(LE,"r")
    MA=zeros(NB,NR+n);
    MA(1,1)=1;

    for i=2:NR
        //CRIAR A MATRIZ A .
        origem=M(i,1);
        destino=M(i,2);
        
        //Parte real matriz incidencia (A).
        MA(origem,i)=-1;
        MA(destino,i)=1;
//        pause
    end
//    pause
endfunction

// ---------------------------------------------
// Funcao responsavel pela criação da Matriz (b).
//----------------------------------------------

function MB=MatrizB(M, LE, vf, vi, buss, b)
    MB=b;
//    if b==[]
//        MB=zeros(NB,1);
//    end
    //disp(vi)
    for i=1:size(buss,"r")
        //CRIAR MATRIZ b a partir das informações obtidas pela matriz Buss.
        aux = buss(i,1); 
        p=M(aux,3); //Potência necessária para subrir a parte real do sistema.
        //Calculo da corrente.
        cu=p/vi(i,1);
        //Parte real atualizada da matriz igualdade b.
        MB(aux)=cu;
    end
    //Inserção do valor de referencia para tensão em cada laço.
    n = size(LE,"r")
    for i=1:n
       MB(i+NB)=vf; 
    end
//    pause
endfunction

// ---------------------------------------------
// Funcao responsavel pela criação da Matriz C.
//----------------------------------------------

function c=MatrizC(M, LE)
    n = size(LE,"r")//Número de laços criados.
    
    //Criação da matriz c.
    c=zeros(NR+n,1)
    for i=1:NB
        if(M(i,7)~=0)
            c(i,1)=1;
        end
    end
endfunction

// ---------------------------------------------
// Estutura principal do Algoritimo para Otimização.
// do fluxo de Corrente.
// -------------------------------------------- 

[LE,buss]=lacoE(M); //Chama a função para criação de laços do circuito.
LE = LE';//Matriz de laço transposta.
A=MatrizA(M,LE);//Chama a função para crianção da matriz A.
c=MatrizC(M,LE);//Chama a função para crianção da matriz c.

dim=1;
A=cat(dim,A,LE);//Faz a concatenação da Matriz A com a Matriz de laços.

DeltaV=e;// Condição inicial de DeltaV
n = size(LE,"r")
vi(1:n,1)=1;// Criação da matriz de tensão dos lacos

xopt=n;
DeltaV2=0;
b=[];
while (DeltaV >= e)
    //Função para atualização e criação da matrizB
    b=MatrizB(M,LE,vf,vi,buss,b);

    //disp("proximo é DeltaV1",DeltaV1)
    DeltaV1=(sum(xopt,"r"))/n;

    //Função de otimização KARMARKAR - objetivo: Minimizar IR no FEEDER
    //Minimize x0 + x12 + x13 + x23 ...Xij
    [xopt,fopt,exitflag,iter,yopt]=karmarkar(A,b,c); 
    yopt.ineqlin
    yopt.lower 
    if(exitflag == 1)
        vi=xopt(NR+1:NR+n,:);
    end
    DeltaV2=(sum(xopt,"r"))/n;
    DeltaV=DeltaV1-DeltaV2;//Declara o novo valor  de Delta
    
    norm(A*xopt-b)
//    pause
end
