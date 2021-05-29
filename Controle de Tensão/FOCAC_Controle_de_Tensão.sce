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
//entradaDeDados=input("Digite o endereço com localização do arquivo de entrada de dados. Obs: seguir instruções no arquivo Guideline. ");

// Arquivo de Entrada com a estrutura da Rede de Distribuição.
//M = fscanfMat(entradaDeDados, "%lg"); // Importa arquivo que apresenta os dados de entrada da rede. 

M = fscanfMat("C:\Users\Higor\Desktop\rede5barras.txt", "%lg");

// Separação de dados a serem usandos ao longo do algoritimo.
// Número de barras da  rede
NBc=M(1,2);
// Salva o número de ramos
NR=M(1,1);
// Salva o número de barras
NB=NBc+M(1,3); 
// Valor de refencia de tensão inicial Real.
vrReal=1.0;
// Valor de referencia para tensão inicial Imag.
vrImag=0;
// Tensão mínima estabelecida para tensão na nova barra de geração
vi=(0.977);

// Inserção de numeros para barra de geração antes iguais a 0.
// Facilitando a criação de um arquivo de entrada para o usuário.
a=1;
for i=2:size(M,"r")
    if M(i,1)==0
        M(i,1)=NBc+a;
        a=a+1;
    end
end

// ---------------------------------------------------------------
// Funcao responsavel pela construção da Matriz incidência A.
// Entrada : Estrutura do arquivo completo e Matriz laços externos.
// Saída : Matriz incidência para restrição de igualdade Aeq
// inequaldade A, informação de colunas para resolução de caso real,
// informação de colunas para resolução de caso completo de rede.
//----------------------------------------------------------------

function [Aeq,qnt_coluna_MA_incidencia,qnt_coluna_MA,aux,A,barraDefeitoFalha] = MatrizA(M,controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao)

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
            // Criação da matriz incidência para restrições de inegualdade.
            A(a,i-1)=1;
            aux=aux+1; // Informação de quantos ramos da rede são de linhas  abertas.
            
            // Criação da matriz referente aos limites da rede, afim de criar uma equação que associa as variáveis inteiras com as não inteiras.
            matrizLimite(a,a)=-5;
            a=a+1; // Atualização de posição da matriz.
        end
    end
//    pause
    // Expansão da matriz de igualdade para inserção das equações que serão responsáveis pela associação de variáveis inteiras e não inteiras.
    Aeq=cat(2,Aeq,zeros(size(Aeq,'r'),aux));
//    pause
    A3=zeros(1,(size(Aeq,'c')-aux)); // Matriz auiliar para criação de restrição inteira
    A3=cat(2,A3,ones(1,aux));
//    A3=cat(2,A3,A4);
//    pause
    //Inserção de restrição para variáveis inteiras
    Aeq=cat(1,Aeq,A3);

    // Comunicação ao usuário para que seja informado qual barra ocorreu o defeito falha.
    desligamento=zeros(1,size(Aeq,"c"));
    barraDefeitoFalha=restricao(M);
    if barraDefeitoFalha==0 then
        continue;
    else
        desligamento(1,barraDefeitoFalha)=1;
    end
    //Inserção da restrição de desligamento para barra.
    Aeq=cat(1,Aeq,desligamento);
    
    // Informação sobre tamanho da matriz incidencia referente a restrições de igualdade.
    qnt_coluna_MA_incidencia=size(Aeq,'c');

    // Ajuste de matriz auxiliar para que seja possível inserir a matriz incidência para restrições de igualdade para parte Imaginária da rede.
    A3=zeros(size(Aeq,'r'),size(Aeq,'c'));
    A3=cat(2,A3,Aeq);
    // Expansão da matriz incidência restrições Real, para que seja inserida restrições para parte imaginária da rede.
    Aeq=cat(2,Aeq,zeros(size(Aeq,'r'),size(Aeq,'c')));
    Aeq=cat(1,Aeq,A3);

    // Expansão da matriz incidência real para inserção das equações de laço criadas.
    Aeq=cat(2,Aeq,zeros(size(Aeq,"r"),size(controleDeTensaoReal,"c")-size(Aeq,"c")));

    // Incremento das esquações de laços criadas para controle de tensão.
    Aeq=cat(1,Aeq,controleDeTensaoReal);
    Aeq=cat(1,Aeq,controleDeTensaoImag);
    
    // Matriz A referente a inequações, inegualdades.
    ANegativo=A*(-1); // Mudança de sinal para criação das restições de inegualdade de forma adequada.
    A=cat(2,A,zeros(size(A,'r'),2*M(1,3))); // Expansão da matriz incidência para restrições de ingualdade
    
    // Adição da matriz de limites à matriz incidência para restrições de ingualdade.
    A=cat(2,A,matrizLimite);
    // Ajuste de tamanho.
    ANegativo=cat(2,ANegativo,zeros(size(ANegativo,'r'),2*M(1,3)));
    // Adição da matriz de limites à matriz negativa de incidência para restrições de ingualdade.
    ANegativo=cat(2,ANegativo,matrizLimite);
    // Junção de todas restrições de inegualdade para associação de variáveis inteiras e não inteiras.
    A=cat(1,A,ANegativo);

    // Inserção para parte Imaginária da rede.
    A2=zeros(size(A,'r'),size(A,'c'));
    A2=cat(2,A2,A);
    
    // Expansão da matriz para inegualdades
    A=cat(2,A,zeros(size(A,'r'),size(A,'c')));
    A=cat(1,A,A2);
    
    // Expansão da matriz de restrições de inegualdades para inserção de equação para definição de limite infeior de tensão de barra.
    A=cat(2,A,zeros(size(A,"r"),size(controleDeTensaoReal,"c")-size(A,"c")));
    // Inserção da restrição de limite inferior para tensão de barra.
    A=cat(1,A,limiteInferiorParaTensao)
    
    // Informação sobre quantidade de colunas da matriz Aeq.
    qnt_coluna_MA=size(Aeq,'c');
endfunction

// --------------------------------------------------------------------
// Funcao responsavel pela contrução da Matriz de carga b.
// Entrada: Estrutura do arquivo completo, informação de colunas para 
// resolução de caso real, informação de colunas para resolução de caso
// completo de rede.
// Saída : Matriz de carga b.
//---------------------------------------------------------------------

function [beq,b] = MatrizB(M,C,controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,barraDefeitoFalha)
    // Definição inicial dos vetores de carga real e imaginário.
    beq=zeros(NB,1); // Real.
    b1eq=zeros(NB,1); // Imaginário.

    // Criação da matriz B com todos valores negativos.
    for i=2:NR+1
        // Verificação da característica da barra, 1=geração, 
        if(M(i,7)==1)
            // Matriz b parte real da geração.
            beq(M(i,1),1)=(M(i,3))//+M((i+1),9));
    
            //Matriz B parte imaginário da geração.
            b1eq(M(i,1),1)=M(i,4);
        elseif(M(i,3)~=0)
            // Matriz b parte real de carga.
            beq(M(i,2),1)=(-1)*(M(i,3));
            // Matriz b parte real de carga.
            b1eq(M(i,2),1)=(-1)*(M(i,4));
        end
    end
    
    // Complemento às cargas da rede, definição resultados esperados para as equações de restrição da rede.
    // Restrição de inteiro para Real.
    beq=cat(1,beq,ones(1,1));
    // Restrição para desligamento de barra Real
    beq=cat(1,beq,zeros(1,1));
    // Junção matriz de carga beq Real com restrições, com Imaginário.
    beq=cat(1,beq,b1eq);
    // Restrição de inteiro para Imaginário.
    beq=cat(1,beq,ones(1,1));
    // Restrição para desligamento de barra Imaginário.
    beq=cat(1,beq,zeros(1,1));
    //Inserção de zeros referente a tensão das equações de laço Real e Imaginário.
    // Inserção de tensão inicial Real de barra.
    beq=cat(1,beq,vrReal*ones(size(controleDeTensaoReal,"r"),1));
    // Inserção de tensão inicial Imaginária de barra.
    beq=cat(1,beq,vrImag*ones(size(controleDeTensaoImag,"r"),1));
    
    // Criação da matriz de restrições para inequações (<=) da rede.
    b=zeros(size(C,'r')-size(limiteInferiorParaTensao,'r'),1);
    // Definição do limite inferior para controle de tensão.
    b=cat(1,b,vi*(-1)*ones(size(limiteInferiorParaTensao,"r"),1));
endfunction

//-------------------------------------------------------------------
// Função responsável pela criação da matriz H.
// Entrada : Estrutura do arquivo completo, informação de colunas para
// resolução de caso real, informação de colunas para resolução de caso
// completo de rede.
// Saída : Matriz simétrica para resistência e reatância H.
//-------------------------------------------------------------------

function [H]=MatrizH(M,qnt_coluna_MA_incidencia,aux,controleDeTensaoReal,controleDeTensaoImag)
    
    AuxG=1;// Auxiliar para criação da matriz inciência de "geração", endereço da coluna.

    //Dimensionamento de tamanho para matriz Q e auxiliares.
    H=zeros(NR,qnt_coluna_MA_incidencia);
    H2=H;
    
    //Criação da matriz simétrica para resistência e reatância de cada barra.
    for i=2:(NR+1)
        //Matriz simétrica para resistência.
        H(i-1,i-1)=M(i,5);
        
        //Matriz simétrica para reatância.
        H2(i-1,i-1)=M(i,6);
    end
    
    //Criação de complemento que será responsável pela resistência e impedância das as barras de geração.
    AuxH=0.0000001*eye(M(1,3),M(1,3));
    
    //Responsável pela criação de uma matriz complementar, responsável para ordenar e possibilitar a concatenação posteriormente (inserir as variáveis responsáveis pela variável de folga).
    H1=zeros(M(1,3),NR);
    
    //Concatenação para que seja possível incrementar as variáveis com baixa resistencia.
    H1=cat(2,H1,AuxH);
    
    // Incremento de restrições para variáveis de restição inteiras Reais.
    H1=cat(2,H1,zeros(M(1,3),aux));
    H=cat(1,H,H1); //Real.
    H3=zeros(aux,size(H,'r'));
    H3=cat(2,H3,zeros(aux,aux));
    H=cat(1,H,H3);
    
    // Incremento de restrições para variáveis de restição inteiras Imaginárias.
    H1=cat(2,H1,zeros(1,size(H2,'c')-size(H1,'c')))
    H2=cat(1,H2,H1); //Real.
    H3=zeros(aux,size(H2,'r'));
    H3=cat(2,H3,zeros(aux,aux));
    H2=cat(1,H2,H3);

    // Junção de pesos para real e imaginário.
    H=cat(2,H,zeros(size(H,'r'),size(H,'r')));
    H2=cat(2,zeros(size(H2,'r'),size(H2,'r')),H2);
    H=cat(1,H,H2);

    // Inserção para pesos das variáveis referentes aos laços da rede.
    H=cat(2,H,zeros(size(H,"r"),size(controleDeTensaoReal,"r")));
    H=cat(2,H,zeros(size(H,"r"),size(controleDeTensaoImag,"r")));
    H=cat(1,H,zeros(size(controleDeTensaoReal,"r"),size(H,"c")));
    H=cat(1,H,zeros(size(controleDeTensaoImag,"r"),size(H,"c")));
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
// Saída : Informação de falha de ramo.
//---------------------------------------------------------------

function [barraDefeitoFalha]=restricao(M)
    while 1>0 do
        // Comunicação ao usuário.
        disp("Houve um Defeito Falha em algum ramo?");
        algumRamoFalhou=input("Digita 1 (um) para sim ou 0 (zero) para não. ");
        if algumRamoFalhou==1
            barraDefeitoFalha=input("Digite o ramo com Defeito Falha: ");
            // Verificação de possibilidade de existir o Defeito Falha informado.
            if(barraDefeitoFalha>NB)
                disp("Você digitou um valor inválido");
            elseif(barraDefeitoFalha~=0)
                // Inserção de valor 1 para iniciar uma restrição capaz de desligar o ramo informado tanto para parte real quanto imaginária.
                break;
            else
                disp("Você digitou um valor inválido");
            end
        else
            barraDefeitoFalha=0;
            break;
        end
    end
endfunction

// --------------------------------------------------------------------
// Funcao responsável pela construçao do laço externo.
// Entrada : Estrutura do arquivo completo.
// Saida : Matriz de lacos externos, informação de quantos laços foram
// criados.
// -------------------------------------------------------------------

function [controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,M]=controleDeTensao(M)
    // Busca a posição das barras onde se encontram os feeders da rede.
    insercaoDeCarga=find(M(:,7)==1)
    
    // Vetor que traz as informações de quais barras estão recebendo potência.
    quaisSaoAsBarrasDeInsercao=M(insercaoDeCarga,2);
    
    // Criação do vetor que será responsável por montar a equação a fim de controlar a tensão Real do novo feeder inserido.
    controleDeTensaoReal=zeros(NR+M(1,3),1);
    
    // Criação do vetor que será responsável por montar a equação a fim de controlar a tensão Imaginária do novo feeder inserido.
    controleDeTensaoImag=zeros(NR+M(1,3),1);
    
    // Definição de onde foi inserida a nova geração.
    novaFonteInserida=M(size(M,"r"),2);
    
    // Posição de todas barras que recebem potência de um feeder.
    posicaoTodasBarrasDeGeracao=find(M(:,7)==1);
    
    // Definição de todas barras que recebem potência de um feeder.
    barrasDeGeracao=M(posicaoTodasBarrasDeGeracao,2);
    
    // Busca pelo feeder mais próximo ao novo feeder inserido
    feedersAntecessores=find(barrasDeGeracao(:,1)<novaFonteInserida);

    // Definição do feeder mais próximo ao novo feeder inserido
    feederMaisProximo=feedersAntecessores(size(feedersAntecessores,"c"));
    
    // Criação da equação de laço para controle de tensão real e imaginário.
    for i=barrasDeGeracao(feederMaisProximo,1):NR+M(1,2)
        if i>=novaFonteInserida
            break
        end
        controleDeTensaoReal(i,1)=M(i+1,5); // Real.
        controleDeTensaoImag(i,1)=M(i+1,6); // Imaginário.
    end
    
    // Transposição das matrizes para controle de tensão real e imaginário.
    controleDeTensaoReal=controleDeTensaoReal'; // Real.
    controleDeTensaoImag=controleDeTensaoImag'; // Imaginário.
    
    // Criação de uma matriz auxiliar para criação dos laços completos.
    controleDeTensaoRealAux=controleDeTensaoReal;
    
    // Posição dos ramos abertos na matriz M
    ramoAberto=find(M(:,8)==1);
    // Informação de quantos ramos são abertos na topologia inicial
    quantosRamosAbertos=size(ramoAberto,"c");

    // Ajuste de dimensões para equações de controle de tensão Real e Imaginário.
    controleDeTensaoReal=cat(2,controleDeTensaoReal,zeros(size(controleDeTensaoReal,"r"),quantosRamosAbertos)); // Real.
    controleDeTensaoImag=cat(2,controleDeTensaoImag,zeros(size(controleDeTensaoImag,"r"),quantosRamosAbertos)); // Imaginário.

    // Criação dos laços completos para controle de tensão.
    controleDeTensaoReal=cat(2,controleDeTensaoReal,controleDeTensaoImag*(-1)); // Real.
    controleDeTensaoImag=cat(2,controleDeTensaoImag,controleDeTensaoRealAux); // Imaginário

    // Expansão do vetor para que possa ser inseridas variáveis para cálculo de tensão.
    controleDeTensaoReal=cat(2,controleDeTensaoReal,zeros(size(controleDeTensaoReal,"r"),size(controleDeTensaoReal,"r")+size(controleDeTensaoImag,"r"))); // Real.
    controleDeTensaoImag=cat(2,controleDeTensaoImag,zeros(size(controleDeTensaoImag,"r"),size(controleDeTensaoReal,"r")+size(controleDeTensaoImag,"r"))); // Imagiário.
    
    // Criação do vetor para restição de tensão mínima na barra de geração.
    limiteInferiorParaTensao=zeros(size(controleDeTensaoReal,"r"),size(controleDeTensaoReal,"c"));
    
    // Inserção de variáveis para calculo de tensão.
    controleDeTensaoReal(1,2*(NR+M(1,3)+quantosRamosAbertos)+1)=1; // Real.
    controleDeTensaoImag(1,2*(NR+M(1,3)+quantosRamosAbertos)+2)=1; // Imaginário.
    
    // Identificação de qual variável deverá controlar a tensão.
    limiteInferiorParaTensao(1,2*(NR+M(1,3)+quantosRamosAbertos)+1)=-1;
endfunction

// --------------------------------------------------
// Estutura principal do Algoritimo para Otimização
// do fluxo de Corrente Alternada com Contingência.
// --------------------------------------------------

// Instrução para criação da matriz responsável pela criação dos laços externos da rede.
[controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,M]=controleDeTensao(M);

// Instrução para criação da matriz incidência A.
[Aeq,qnt_coluna_MA_incidencia,qnt_coluna_MA,aux,A,barraDefeitoFalha]=MatrizA(M,controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao);

// Instrução para criação da matriz incidência Q.
[H]=MatrizH(M,qnt_coluna_MA_incidencia,aux,controleDeTensaoReal,controleDeTensaoImag);

// Instrução para criação da matriz incidência P.
[f]=MatrizF(qnt_coluna_MA);

// Instrução para criação da matriz incidência B.
[beq,b]=MatrizB(M,A,controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,barraDefeitoFalha);

// Definição das variáveis inteiras/binárias da rede.
// Encontrar as variáveis inteiras/binárias reais da rede.
intconReal=find(Aeq(NB+1,:)==1);
// Encontrar as variáveis inteiras/binárias imaginárias da rede.
intconImaginario=find(Aeq(2*(NB+1)+1,:)==1);
// Junção das variáveis inteiras/binárias reais e imaginárias da rede.
// As variáveis reais e imaginárias da rede serão as mesmas, apenas deslocadas umas das outras.
intcon=cat(2,intconReal,intconImaginario);

// Função de otimização FOT_INTQUADPROG - objetivo: Minimização de perdas na rede com opções de manobra.
[xopt,fopt,exitflag,output]=fot_intquadprog(H,f,intcon,A,b,Aeq,beq)

// Comunicação ao usuário sobre os resultados obtidos com a função de otimização.
if exitflag == 0 then
    disp("Solução Ótima Encontrada!");
    disp(fopt,"O valor ótimo encontrado para a função objetivo.");
elseif exitflag == 1 then
    disp("Solução não encontrada.")
else
    disp("Erro encontrado.")
end
