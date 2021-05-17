//----------------------------------------------
// Algorítimo OCF responsável pela otimização na
// corrente injetada na rede distribuição.
///
// @date 01/07/2020
// @author Higor de Paula Kolecha
// @author Adolfo Blengini Neto
// @author Marcius Fabius Henriques de Carvalho.
// @version 1.0
//----------------------------------------------
clear;
clc;

// Arquivo de Entrada com a estrutura da Rede de Distribuição.
M = fscanfMat("C:\Users\Higor\Desktop\Iniciação científica\IC 2020 marcius scilab\Programação\11 barras v2\inputnelson.txt", "%lg"); // Importa arquivo que apresenta os dados de entrada da rede. 

// Salva o número de ramos
NR=M(1,1); 

// Salva o número de barras
NB=M(1,2); 

//valor de refencia de tensão inicial real
vf=1.02;

//valor de refencia de tensão inicial imaginario
vfi=0;

//valor de convergencia
e=0.0001;

// -----------------------------------------------------
// Funcao responsável pela construçao do laço externo 
// Entrada : Estrutura do arquivo completo 
// Saida : Matriz de lacos externos
// -----------------------------------------------------

function [quaislacos,barralaco] = lacoE (M)
    //Entrada
    M(1,1)=1; // Inicia os contadores de linhas e colunas
    M(1,2)=1; // Inicia os contadores de linhas e colunas
    
    O=M(:,1);  // Criação do vetor origem.
    D=M(:,2);  // Criação do vetor destino.

    quaislacos=[] // Criação do caminho que a corrente faz para poder criar um laço.
    quaislacosimaginario=[]
    barralaco=[]
    barralacoimaginario=[]
    for i=NB:-1:1
        buss=i;     
        if (M(i,3)~=0) //Verificação da existencia de carga. Caso afirmativo, passa para a criação dos laços
            
            laco=zeros(NB,1); // Criação do caminho que a corrente faz para poder criar um laço.
            
            barralaco=cat(1,barralaco,i);// Faz a concatenação da matriz "barralaco" a fim de serem adicionadas apenas laços existentes
            
            //Busca amontante até feeder, ou seja, barra 1
            while (buss~=1)
                o=find(D(:,1)==buss); //Mostra a posição onde se encontra o próximo ramo para formação de laços
                laco(o,1)=M(o,5); //Monta o vetor caminhos, mostrando os laços feitos
                buss=O(o); //Recebe o valor que se encontra na posição acima destacada dentro do vetor origem
            end
            buss=i;
            quaislacos = cat(2,quaislacos,laco);
        end

        if(M(i,4)~=0)
            lacoimag=zeros(NB,1);
            barralacoimaginario=cat(1,barralacoimaginario,i);
            while (buss~=1)
                o=find(D(:,1)==buss);
                lacoimag(o,1)=M(o,6);
                buss=O(o);
            end
            quaislacosimaginario=cat(2,quaislacosimaginario,lacoimag);
         end
    end
    quaislacos2=[];
    quaislacosimaginario2=[];
    
    barralaco=cat(1,barralaco,barralacoimaginario);
// -----------------------------------------------------------------------------
// Multiplicação das matrizes por -1.
// -----------------------------------------------------------------------------
    for i=1:(size(quaislacos,"r"))
        for j=1:(size(quaislacos,"c"))
            if(quaislacos(i,j)~=0)
                quaislacos2(i,j)=quaislacos(i,j)*(-1);
            end
            if(quaislacosimaginario(i,j)~=0)
                quaislacosimaginario2(i,j)=quaislacosimaginario(i,j)*(-1);
            end
        end
    end
    
    quaislacos=cat(2,quaislacos,quaislacos2);
    
    quaislacosimaginario=cat(2,quaislacosimaginario2,quaislacosimaginario);
    
    quaislacos=cat(1,quaislacos,quaislacosimaginario);
    
    l=size(quaislacos,"c")// Mede o tamanho da matriz quaislaacos
    complemento=eye(l,l);
    quaislacos=cat(1,quaislacos,complemento);// Faz a concatenação da matriz quaislacos com a matriz complemento
    
endfunction //Fim da função de criação de laços

// ---------------------------------------------
// Funcao responsavel pela criação da Matriz A
//----------------------------------------------

function MA = MatrizA(M,LE)
    
    n = size(LE,"r")

    MA=zeros(2*NB,2*NR+n);
    
    MA(1,1)=1;MA(1+NR,1+NR)=1;
    for i=2:NR
        //CRIAR A MATRIZ A .
        origem=M(i,1);
        destino=M(i,2);
        
        //Parte real matriz incidencia (A).
        MA(origem,i)=-1;
        MA(destino,i)=1;
        
        //Parte imaginária da matriz incidencia (A)
        MA(NB+origem,i+NB)=-1;
        MA(NB+destino,i+NB)=1;
    end

endfunction

// ---------------------------------------------
// Funcao responsavel pela criação da Matriz (b).
//----------------------------------------------

function MB = MatrizB(M,LE,vf,vi,buss,b)

    MB=b;
    qntlaco=((size(buss,"r"))/2);
    for i=1:qntlaco
        //CRIAR MATRIZ b a partir das informações obtidas pela matriz Buss.
        aux = buss(i,1); 
        p=M(aux,3); //Potência necessária para subrir a parte real do sistema.
        
        //Calculo da corrente.
        cu=p/vi(i,1);
        
        //Parte real atualizada da matriz igualdade b.
        MB(aux)=cu;
        
        //imaginario
        aux2=buss(i+qntlaco,1);
        p=M(aux2,4);
        
        cu=p/vi((i+qntlaco),1);
        
        MB(aux2+NB)=(cu*(1));
    end
    
    //Inserção do valor de referencia para tensão em cada laço.
    for i=1:qntlaco
       MB((2*NR+i),1)=vf;
       MB((2*NR+qntlaco+i),1)=vfi;
    end
endfunction

// ---------------------------------------------
// Funcao responsavel pela criação da Matriz C.
//----------------------------------------------

function c=MatrizC(A,LE)
    
    n = size(A,"c")//Número de laços criados.

    //Criação da matriz c.
    c=zeros(n,1)
    
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


A=cat(1,A,LE);//Faz a concatenação da Matriz A com a Matriz de laços.

c=MatrizC(A,LE);//Chama a função para crianção da matriz c.

n = size(LE,"r")
vi(1:n,1)=1;// Criação da matriz de tensão dos lacos

b=[];

//Função para atualização e criação da matrizB
b=MatrizB(M,LE,vf,vi,buss,b);
    
//Função de otimização KARMARKAR - objetivo: Minimizar IR no FEEDER
//Minimize x0 + x12 + x13 + x23 ...Xij
[xopt,fopt,exitflag,iter,yopt]=karmarkar(A,b,c); 
yopt.ineqlin
yopt.lower 

norm(A*xopt-b)
