//-----------------------------------------------
// Algorítimo FCCC responsável pela otimização na
// corrente injetada na rede distribuição.
// Projeto de Pesquisa FAPESP.
// Projeto número: #2019/24128-2.
// @date 01/07/2020.
// @author Higor de Paula Kolecha.
// @author Adolfo Blengini Neto.
// @author Marcius Fabius Henriques de Carvalho.
// @version 1.0
//-----------------------------------------------

// Responsável pela limpeza toda memória
clear;
// Responsável pela limpeza de tela
clc;

// Solicitação ao usuário do endereço para obtenção do arquivo de entrada.
entradaDeDados=input("Digite o endereço com localização do arquivo de entrada de dados. Obs: seguir instruções no arquivo Guideline. ");

// Arquivo de Entrada com a estrutura da Rede de Distribuição.
M = fscanfMat(entradaDeDados, "%lg"); // Importa arquivo que apresenta os dados de entrada da rede. 

// Inicio de contagem de timer para tempo de resolução da rede.
tic();

// Salva o número de ramos
NR=M(1,1); 
// Salva o número de barras
NB=M(1,2); 
// Valor de refencia de tensão inicial
vf=1.0;
// Valor de convergencia
e=0.0001;

// -----------------------------------------------------
// Funcao responsável pela construçao do laço externo 
// Entrada : Estrutura do arquivo completo 
// Saida : Matriz de lacos externos
// -----------------------------------------------------

function [quaislacos, barralaco]=lacoE(M)
    // Entrada de dados.
    M(1,1)=1; // Inicia os contadores de linhas e colunas
    M(1,2)=1; // Inicia os contadores de linhas e colunas
    
    O=M(:,1); // Criação do vetor origem.
    D=M(:,2); // Criação do vetor destino.
    
    quaislacos=[] // Criação do caminho que a corrente faz para poder criar um laço.
    
    barralaco=[];
    for i=NR:-1:1
        buss=i;     
        if (M(i,3)~=0) // Verificação da existencia de carga. Caso afirmativo, passa para a criação dos laços.
            laco=zeros(NR,1); // Criação do caminho que a corrente faz para poder criar um laço.
            barralaco=cat(1,barralaco,i); // Faz a concatenação da matriz "barralaco" a fim de serem adicionadas apenas laços existentes.

            // Busca amontante até feeder, ou seja, barra 1
            while (buss~=1)
                o=find(D(:,1)==buss); // Mostra a posição onde se encontra o próximo ramo para formação de laços.
                laco(o(1,1),1)=M(o(1,1),5); // Monta o vetor caminhos, mostrando os laços feitos.
                buss=O(o(1,1)); // Recebe o valor que se encontra na posição acima destacada dentro do vetor origem.
            end
            // Criação da matriz final para criação de laços.
            quaislacos = cat(2,quaislacos,laco);
        end
    end
    // Informação sobre o tamanho da matriz criada referente a quantidade de colunas.
    l = size(quaislacos,"c")
    
    // Criação da matriz complemento para definição de variáveis.
    complemento = eye(l,l);
    
    // Junção da matriz complemento na matriz original de laços criados.
    quaislacos = cat(1,quaislacos, complemento); 
    
    // Transposição da matriz quaislacos.
    quaislacos = quaislacos';
endfunction //Fim da função de criação de laços

// ---------------------------------------------------------------
// Funcao responsavel pela construção da Matriz incidência A.
// Entrada : Estrutura do arquivo completo e Matriz laços externos.
// Saída : Matriz incidência.
//----------------------------------------------------------------

function MA=MatrizA(M, LE)
    // Determinação da quantidade de laços criados, somando-se laço referente a parte real e imaginária da rede.
    n = size(LE,"r");
    // Dimensionamento da matriz incidência MA.
    MA=zeros(NB,NR+n);
    // Determinação da barra de geração, neste caso a barra 1 está recebendo a carga.
    MA(1,1)=1;
    
    // Estrutura de repetição para criação da matriz incidência MA.
    for i=2:NR
        // CRIAR A MATRIZ A .
        origem=M(i,1);
        destino=M(i,2);
        
        // Parte real matriz incidencia (A).
        MA(origem,i)=-1;
        MA(destino,i)=1;
    end
    
    // Junção da Matriz Incidência com matriz de Laços.
    MA=cat(1,MA,LE);
endfunction // Fim da função criação de matriz incidência.

// --------------------------------------------------------------------
// Funcao responsavel pela contrução da Matriz de carga b.
// Entrada: Estrutura do arquivo completo, Matriz laços externos,
// Valor de referência de tensão inicial, Tensão obtida pelo algoritmo,
// informação de quais barras criaram laços e matriz de carga antiga.
// Saída : Matriz de carga atualizada b.
//---------------------------------------------------------------------

function MB=MatrizB(M, LE, vf, vi, buss, b)
    // Definição inicial da matriz carga B
    MB=b;
    // Atualização dos valores de carga para convergencia.
    for i=1:size(buss,"r")
        // CRIAR MATRIZ b a partir das informações obtidas pela matriz Buss.
        aux = buss(i,1); 
        p=M(aux,3); // Potência necessária para subrir a parte real do sistema.
        // Calculo da corrente.
        cu=p/vi(i,1);
        // Parte real atualizada da matriz igualdade b.
        MB(aux)=cu;
    end
    
    // Inserção do valor de referencia para tensão em cada laço.
    n = size(LE,"r")
    for i=1:n
       MB(i+NB)=vf; 
    end
endfunction // Fim da criação da matriz de carga.

// ----------------------------------------------------------------
// Funcao responsavel pela criação da Matriz c.
// Entrada : Estrutura do arquivo completo e Matriz laços externos.
// Saída : Matriz peso c.
//-----------------------------------------------------------------

function c=MatrizC(M, LE)
    n = size(LE,"r"); // Número de laços criados.
    // Dimensionamento da matriz peso c.
    c=zeros(NR+n,1)
    // Estrutura de repetição para criação da matriz c.
    for i=1:NB
        if(M(i,7)~=0)
            c(i,1)=1;
        end
    end
endfunction // Fim da função criação de pesos da função.

// ------------------------------------------------
// Estutura principal do Algoritimo para Otimização.
// do fluxo de Corrente Contínua.
// ------------------------------------------------

[LE,buss]=lacoE(M); // Chama a função para criação de laços do circuito.
A=MatrizA(M,LE); // Chama a função para crianção da matriz A.
c=MatrizC(M,LE); // Chama a função para crianção da matriz c.

DeltaV=e; // Condição inicial de DeltaV.
n = size(LE,"r"); // Determinação da quantiade de laços criados.
vi(1:n,1)=1; // Criação da matriz de tensão dos lacos.

xopt=n; // Determinação inicial de qualquer valor.
// Determinação de valores iniciais.
DeltaV2=0;
b=[];

// Definição dos limites superiores (ub) e inferiores (lb) dos ramos da rede.
lb=(-5)*ones(size(A,"c"),1);
ub=5*ones(size(A,"c"),1);

// Definição dos valores iniciais para as variáveis da rede.
x0=[];

// Estrutura de repetição a fim de chegar em convergencia
while (DeltaV >= e)
    // Função para atualização e criação da matrizB
    b=MatrizB(M,LE,vf,vi,buss,b);

    DeltaV1=(sum(xopt,"r"))/n; // Calculo para DeltaV1.

    // Função de otimização KARMARKAR - objetivo: Minimizar IR no FEEDER
    // Minimize x0 + x12 + x13 + x23 ...Xij
    [xopt,fopt,exitflag,iter,yopt]=karmarkar(A,b,c,x0,[],[],[],[],[],[],lb,ub); 
    yopt.ineqlin
    yopt.lower 
    // Atualização dos valores de tensão caso a função ache uma otimização.
    if(exitflag == 1)
        vi=xopt(NR+1:NR+n,:);
    end
    DeltaV2=(sum(xopt,"r"))/n; // Cálculo para valor de média de resultados para DeltaV2.
    DeltaV=DeltaV1-DeltaV2; // Declara o novo valor  de Delta

    norm(A*xopt-b)
end

// Término de contagem de timer para tempo de resolução da rede.
toc();

// Informação ao usuário da resolução da rede.
if xopt~=[]
    disp("Solução Ótima Encontrada!");
    disp(fopt,"O valor ótimo encontrado para a função objetivo.");
    disp(iter,"Foram necessárias esta quantidade de iterações para convergir.");
    disp("O primeiro valor se refere às iterações e o segundo às restrições desativadas para resolução da rede.");
    disp(ans,"CPU time (s).");
else
    disp("Solução não encontrada.");
end
