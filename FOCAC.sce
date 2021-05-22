//----------------------------------------------
// Algorítimo OCF responsável pela otimização na
// corrente injetada na rede distribuição.
//Projeto FAPESP
//nº...
///
// @date 01/07/2020
// @author Higor de Paula Kolecha
// @author Adolfo Blengini Neto
// @author Marcius Fabius Henriques de Carvalho.
// @version 1.2
//----------------------------------------------

// Responsável pela limpeza toda memória.
clear;
// Responsável pela limpeza de tela.
clc;

//Solicitação ao usuário do endereço para obtenção do arquivo de entrada.
entradaDeDados=input("Digite o endereço com localização do arquivo de entrada de dados. Obs: seguir instruções no arquivo Guideline. ");

// Arquivo de Entrada com a estrutura da Rede de Distribuição.
M = fscanfMat(entradaDeDados, "%lg"); // Importa arquivo que apresenta os dados de entrada da rede. 

// Salva o número de ramos
NR=M(1,1);
// Salva o número de barras
NB=M(1,2)+M(1,3); 
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

// ---------------------------------------------------------------
// Funcao responsavel pela construção da Matriz incidência A.
// Entrada : Estrutura do arquivo completo e Matriz laços externos.
// Saída : Matriz incidência para restrição de igualdade Aeq
// inequaldade A, informação de colunas para resolução de caso real,
// informação de colunas para resolução de caso completo de rede.
//----------------------------------------------------------------

function [Aeq,qnt_coluna_MA_incidencia,qnt_coluna_MA,aux,A] = MatrizA(M,LE)

    AuxC=1; // Auxiliar para criação da matriz incidência de "carga", endereço da coluna.
    AuxG=1; // Auxiliar para criação da matriz incidência de "geração", endereço da coluna.
    aux=0;
    
    MAg=zeros(M(1,2),M(1,3)); // Criação de vaiáveis novas para descartar restrições de desigualdade.
        
    for i=2:(NR+1)
        // Orientação.
        origem=M(i,1);
        destino=M(i,2);
    
        // Parte real matriz incidencia (A).
        Aeq(origem,AuxC)=1;
        Aeq(destino,AuxC)=-1;
        // Atualização da variável auxiliar
        AuxC=AuxC+1;
        
        //Complemento da matriz incidencia (A) com as barras de geração.
        if(M(i,7)==1)
            MAg(origem,AuxG)=1;
            AuxG=AuxG+1;
        end
    end
    
    //Concatenação da Matriz indicencia com variáveis de folga.
    Aeq=cat(2,Aeq,MAg);

    //Criação da matriz que será responsável por receber quais ramos devem estar normalmente abertas ou fechadas.
    A=[];    
    
    //acrescentado para inteiro
    matrizLimite=[];
    //-------------------------
    
    a=1;
    
    for i=2:(NR+1)
        if M(i,8)==1
            A(a,i-1)=1;
            aux=aux+1;
            //acrescentado para inteiro
            matrizLimite(a,a)=-5;
            a=a+1;
            //-------------------------
        end
    end

    //acrescentado para inteiro
    //real
    Aeq=cat(2,Aeq,zeros(size(Aeq,'r'),aux));
    A3=zeros(1,(size(Aeq,'c')-aux));
    A4=ones(1,aux);
    A3=cat(2,A3,A4);
    Aeq=cat(1,Aeq,A3);

    //Teste desligamento barra real
    desligamento=zeros(1,size(Aeq,"c"));
    desligamento(1,9)=1;
    Aeq=cat(1,Aeq,desligamento);
    //-------------------------------

    qnt_coluna_MA_incidencia=size(Aeq,'c');
    
    //imag
    A3=zeros(size(Aeq,'r'),size(Aeq,'c'));
    A3=cat(2,A3,Aeq);
    Aeq=cat(2,Aeq,zeros(size(Aeq,'r'),size(Aeq,'c')));
    Aeq=cat(1,Aeq,A3);
    //----

    //equações de tensão
    Aeq=cat(2,Aeq,zeros(size(Aeq,"r"),size(LE,"c")-size(Aeq,"c")));
    //--------------

    //Adição das equações de laço
    Aeq=cat(1,Aeq,LE);
    //--------------------------

    qnt_coluna_MA=size(Aeq,'c');
    
    //Matriz C, referente a inequações, inegualdades.
    ANegativo=A*(-1);
    A=cat(2,A,zeros(size(A,'r'),M(1,3)));
    A=cat(2,A,matrizLimite);
    ANegativo=cat(2,ANegativo,zeros(size(ANegativo,'r'),M(1,3)));
    ANegativo=cat(2,ANegativo,matrizLimite);
    A=cat(1,A,ANegativo);

    //imag
    A2=zeros(size(A,'r'),size(A,'c'));
    A2=cat(2,A2,A);
    //    
    A=cat(2,A,zeros(size(A,'r'),size(A,'c')));
    A=cat(1,A,A2);
    
    //equações de tensão
    A=cat(2,A,zeros(size(A,"r"),size(LE,"c")-size(A,"c")));
    //------------------

    qnt_coluna_MA=size(Aeq,'c');

endfunction

//-----------------------
//Criação matriz b
//-----------------------
function [beq,b] = MatrizB(M,quantosLacos,C,LE)
    //Criação da matriz B com todos valores negativos.
    for i=1:NB
        //Matriz b real carga.
        beq(i)=(-1)*(M(i+1,3))//+M((i+1),9));

        //Matriz B imaginário carga.
        b1(i)=(-1)*M(i+1,4);
    end

    //Atualização da matriz B para que os valores referentes a geração estejam positivos e os referentes à carga estejam negativos.
    for i=1:NB
        if(M(i+1,7)==1)
            origem=M(i+1,1);
            //Correção parte Real.
            beq(origem)=beq(origem)*(-1);

            //Correção parte imaginária.
            b1(origem)=b1(origem)*(-1);
        end
    end
    
    beq(NB+1,1)=1;
    
    //Teste desligamento barra 5 real
    beq=cat(1,beq,zeros(1,1));
    //-------------------------------

    beq=cat(1,beq,b1);
    
    beq(2*(NB+1)+1,1)=1;
    
    //Teste desligamento barra 5 real
    beq=cat(1,beq,zeros(1,1));
    //-------------------------------
//pause
    //Inserção de zeros referente a tensão das equações de laço real
    beq=cat(1,beq,zeros(size(LE,"r"),1));
    //--------------------------------------------------------------

    b=zeros(size(C,'r'),1);
endfunction

function [Q]=MatrizH(M,qnt_coluna_MA_incidencia,aux,LE)
//-----------------------
//Criação matriz Q.
//-----------------------
    
    AuxG=1;// Auxiliar para criação da matriz inciência de "geração", endereço da coluna.

    //Dimensionamento de tamanho para matriz Q.
    Q=zeros(NR,qnt_coluna_MA_incidencia);
    Q2=zeros(NR,qnt_coluna_MA_incidencia);
    Q3=zeros(NR,qnt_coluna_MA_incidencia);
    
    //Criação da matriz simétrica para resistência e reatância de cada barra.
    for i=2:(NR+1)
        //Matriz simétrica para resistência.
        Q(i-1,i-1)=M(i,5);
        
        //Matriz simétrica para reatância.
        Q2(i-1,i-1)=M(i,6);
    end
    
    //Criação de complemento que será responsável pela resistência e impedância das as barras de geração.
    AuxQ=0.0000001*eye(M(1,3),M(1,3));
    
    //Responsável pela criação de uma matriz complementar, responsável para ordenar e possibilitar a concatenação posteriormente (inserir as variáveis responsáveis pela variável de folga).
    Q1=zeros(M(1,3),NR);
    
    //Concatenação para que seja possível incrementar as variáveis com baixa resistencia.
    Q1=cat(2,Q1,AuxQ);
    
    //acrescentado para inteiro
    Q1=cat(2,Q1,zeros(M(1,3),aux));
    Q=cat(1,Q,Q1); //Real.
    Q6=zeros(aux,size(Q,'r'));
    Q6=cat(2,Q6,eye(aux,aux));
    Q=cat(1,Q,Q6);
    
    //imag
    Q1=cat(2,Q1,zeros(1,size(Q2,'c')-size(Q1,'c')))
    Q2=cat(1,Q2,Q1); //Real.
    Q6=zeros(aux,size(Q2,'r'));
    Q6=cat(2,Q6,eye(aux,aux));
    Q2=cat(1,Q2,Q6);

    //junção real e imag
    Q=cat(2,Q,zeros(size(Q,'r'),size(Q,'r')));
    Q2=cat(2,zeros(size(Q2,'r'),size(Q2,'r')),Q2);
    Q=cat(1,Q,Q2);
    //------------------

    //variaveis a mais
    Q=cat(2,Q,zeros(size(Q,"r"),size(LE,"r")));
    Q=cat(1,Q,zeros(size(LE,"r"),size(Q,"c")));
    //----------------

endfunction

//----------------------------------------------------
// Função responsável pela criação matriz da f.
// Entrada : Valor de colunas da matriz incidência A.
// Saída : Matriz f.
//----------------------------------------------------
function [f]=MatrizF(qnt_coluna_MA)
    // Matriz de coeficientes dos termos lineares no problema quadrático.
    for i=1:qnt_coluna_MA
        f(i,1)=0;
    end
endfunction

//----------------------------------------------------------------
// Função responável por inserção de restrições para abertura e 
// fechamento de ramos.
// Entrada : Matriz incidência C e Estrutura do arquivo completo.
// Saída : Matriz incidência C  com restrições.
//---------------------------------------------------------------

function [C]=restricao(C,M)
    while 1>0 do
        // Comunicação ao usuário.
        restricao=input("Digite o ramo com Defeito Falha: ");
        // Verificação de possibilidade de existir o Defeito Falha informado.
        if(restricao>NB)
            disp("Você digitou um valor inválido");
        elseif(restricao~=0)
            // Inserção de valor 1 para iniciar uma restrição capaz de desligar o ramo informado tanto para parte real quanto imaginária.
            C(2*NB+M(1,3)+restricao,M(1,3)+restricao)=1;
            C(2*NB+NR+(M(1,3)+restricao,NR+restricao+M(1,3)))=1;
            disp(2*NB+M(1,3)+NR+restricao,NR+restricao+M(1,3));
        elseif(restricao==0)
            disp("Você digitou um valor inválido");
        end
    end
endfunction

// --------------------------------------------------------------------
// Funcao responsável pela construçao do laço externo.
// Entrada : Estrutura do arquivo completo.
// Saida : Matriz de lacos externos, informação de quantos laços foram
// criados.
// -------------------------------------------------------------------

function [LE,quantosLacos]=lacosExternos(M)
    mAux=M(2:NR+1,1:size(M,'c'));
    
    O=mAux(:,1);// Criação do vetor origem.
    D=mAux(:,2);// Criação do vetor destino.
    
    quaisLacos=[]; // Criação do caminho que a corrente faz para poder criar um laço para a parte real da rede.
    barraLaco=[]; //Quais barras tiveram seus laços criados para a parte real da rede.
    
    quaisLacosImag=[]; // Criação do caminho que a corrente faz para poder criar um laço para a parte imaginária da rede.
    barraLacoImag=[]; //Quais barras tiveram seus laços criados para a parte imaginária da rede.
    
    geradores=find(mAux(:,7)==1); //Defini a posição na matriz 'teste', onde estão as barras geradoras
    geradores=mAux(geradores,1); // Cria uma matriz com todas barras geradoras da rede.
    
    for i=NR:-1:1
        buss=mAux(i,2); //Defini a barra inicial que será verificada a existência de laços.
    
        //Laços para Potência, ou seja, parte real da rede.
        if (mAux(i,3)~=0) //Verificação da existencia de carga. Caso afirmativo, passa para a criação dos laços
            laco=zeros(NR+M(1,3),1); // Criação do caminho que a corrente faz para poder criar um laço.
            barraLaco=cat(1,barraLaco,i);// Faz a concatenação da matriz "barralaco" a fim de serem adicionadas apenas as barras com laços existentes
//    pause
            //Busca amontante até feeder, ou seja, barra de geração
            while (1>0)
//                pause
                origem=find(D(:,1)==buss); //Mostra a posição onde se encontra o próximo ramo para formação de laços.
                laco(origem(1,1),1)=mAux(origem(1,1),5); //Monta o vetor caminhos, mostrando os laços feitos.
                buss=O(max(origem,"c")); //Recebe o valor que se encontra na posição acima destacada dentro do vetor origem.
                parada=find(O(:,1)==buss(1,1)) //Procura onde está a origem do buss, retornando a posição da matriz O.
                parou=mAux(parada,1) //Atribui o valor referente a posição da matriz O encontrada na linha superior.
                condicaoDeParada=find(geradores(:,:)==parou) //procura na matriz geradores, se já chegou em algum deles
//                pause
                if(condicaoDeParada~=[])//Caso o valor da condicaoDeParada seja diferente de vazio, quer dizer que ainda não estamos em um gerador, desta forma, a função while deverá continuar. Caso contrário, deverá parar.
                    break
                end
            end
            buss=mAux(i,2); //Retorna o valor inicial referente a barra a ser utilizada neste momento, para que possa ser criados outros laços futuramente.
            quaisLacos=cat(2,quaisLacos,laco); //Concatenação horizontal, inserindo todos laços 'reais' criados em apenas uma matriz.
        end
        
        //Laços para Reatancia, ou seja, parte imaginária da rede.
        if (mAux(i,4)~=0) //Verificação da existencia de carga. Caso afirmativo, passa para a criação dos laços
            lacoImag=zeros(NR+M(1,3),1); // Criação do caminho que a corrente faz para poder criar um laço.
            barraLacoImag=cat(1,barraLacoImag,i);// Faz a concatenação da matriz "barralaco" a fim de serem adicionadas apenas as barras com laços existentes
    
            //Busca amontante até feeder, ou seja, barra de geração
            while (1>0)
                origem=find(D(:,1)==buss); //Mostra a posição onde se encontra o próximo ramo para formação de laços.
                lacoImag(origem(1,1),1)=mAux(origem(1,1),6); //Monta o vetor caminhos, mostrando os laços feitos.
                buss=O(origem(1,1),1); //Recebe o valor que se encontra na posição acima destacada dentro do vetor origem.
                parada=find(O(:,1)==buss(1,1)) //Procura onde está a origem do buss, retornando a posição da matriz O.
                parou=mAux(parada,1) //Atribui o valor referente a posição da matriz O encontrada na linha superior.
                condicaoDeParada=find(geradores(:,:)==parou) //procura na matriz geradores, se já chegou em algum deles

                if(condicaoDeParada~=[])//Caso o valor da condicaoDeParada seja diferente de vazio, quer dizer que ainda não estamos em um gerador, desta forma, a função while deverá continuar. Caso contrário, deverá parar.
                    break
                end
            end
            buss=mAux(i,2); //Retorna o valor inicial referente a barra a ser utilizada neste momento, para que possa ser criados outros laços futuramente.
            quaisLacosImag=cat(2,quaisLacosImag,lacoImag);//Concatenação horizontal, inserindo todos laços 'imaginários' criados em apenas uma matriz.
        end
    end
    aux=0;
    for i=2:(NR+1)
        if M(i,8)==1
            aux=aux+1;
        end
    end
    
    quaisLacos2=quaisLacos;

    barraLaco=cat(1,barraLaco,barraLacoImag);

    quaisLacos=cat(2,quaisLacos,quaisLacosImag);
    
    //inserção das variáveis vazias para criação da matriz LE
    quaisLacos=cat(1,quaisLacos,zeros(aux,size(quaisLacos,'c')));
    //-------------------------------------------------------
    
    quaisLacosImag=cat(2,quaisLacosImag*(-1),quaisLacos2);
    quaisLacosImag=cat(1,quaisLacosImag,zeros(aux,size(quaisLacosImag,'c')));

    
    quaisLacos=cat(1,quaisLacos,quaisLacosImag);
    
    quaisLacos=quaisLacos'
    
    quaisLacos=cat(2,quaisLacos,eye(size(quaisLacos,'r'),size(quaisLacos,'r')));
    
    LE=quaisLacos;
    quantosLacos=size(LE,'r');
endfunction

// ------------------------------------------------
// Estutura principal do Algoritimo para Otimização.
// do fluxo de Corrente Alternada com Contingência.
// ------------------------------------------------

//Instrução para criação da matriz responsável pela criação dos laços externos da rede.
[LE,quantosLacos]=lacosExternos(M);

//Instrução para criação da matriz incidência A.
[Aeq,qnt_coluna_MA_incidencia,qnt_coluna_MA,aux,A]=MatrizA(M,LE);

//Instrução para criação da matriz incidência Q.
[H]=MatrizH(M,qnt_coluna_MA_incidencia,aux);

//Instrução para criação da matriz incidência P.
[f]=MatrizF(qnt_coluna_MA);

//[C]=restricao(C,M);

//Instrução para criação da matriz incidência B.
[beq,b]=MatrizB(M,quantosLacos,A);

intcon=[9,10,19,20];
//intcon=[9,10];

//intcon=[28,29,30,31,32,33,34,35];
//intcon=[28,29,30,31,32,33,34,35,63,64,65,66,67,68,69,70];

[xopt,fopt,exitflag,output]=fot_intquadprog(H,f,intcon,A,b,Aeq,beq)

if exitflag == 0 then
    disp("Solução Ótima Encontrada!")
    disp(fopt,"O valor ótimo encontrado para a função objetivo.")
elseif exitflag == 1 then
    disp("Solução não encontrada.")
else
    disp("Erro encontrado.")
end
