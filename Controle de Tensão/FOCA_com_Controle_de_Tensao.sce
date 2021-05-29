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
M = fscanfMat("C:\Users\Higor\Desktop\rede400barras.txt", "%lg"); // Importa arquivo que apresenta os dados de entrada da rede. 

//temNovaGeracao=input("Fora inserida alguma geração a mais na rede? Digite 1 para sim, 0 para não. ");
//
//qualBarra=input("Qual barra está recebendo a nova fonte? ");
//geracaoAtiva=input("Quanto ela gera de potência ativa? ");
//geracaoReativa=input("Quanto ela gera de potência reativa? ");
//resitenciaRamo=input("Qual a resistencia a ser considerada. obs:~=0. ");
//reatanciaRamo=input("Qual a reatancia a ser considerada. obs:~=0. ");
//novaLinha=[0 qualBarra geracaoAtiva geracaoReativa resitenciaRamo reatanciaRamo 1 0 0];
//M=cat(1,M,novaLinha);
//M(1,1)=M(1,1)+1;
//pause

//"C:\Users\Higor\Desktop\Iniciação científica\2021\Rede 34 barras\rede34barras.txt"
//"C:\Users\Higor\Desktop\rede16barras.txt"
//"C:\Users\Higor\Desktop\rede5barras.txt"
// Arquivo de Entrada com a estrutura da Rede de Distribuição.
//M = fscanfMat(entradaDeDados, "%lg"); // Importa arquivo que apresenta os dados de entrada da rede.

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

fonteAMais=[];

a=1;
for i=2:size(M,"r")
    if M(i,1)==0
        M(i,1)=NBc+a;
        a=a+1;
    end
end

// --------------------------------------------------------------------
// Funcao responsável pela construçao do laço externo.
// Entrada : Estrutura do arquivo completo.
// Saida : Matriz de lacos externos, informação de quantos laços foram
// criados.
// -------------------------------------------------------------------

function [controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,M,fonteAMais]=controleDeTensao(M)
//    fonteAMais=input("Gostaria de inserir alguma fonte de geração a mais? 1=sim,2=não. ");
//    if fonteAMais==1 then
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
        
        // Criação dos laços completos para controle de tensão.
        controleDeTensaoReal=cat(2,controleDeTensaoReal,controleDeTensaoImag*(-1)); // Real.
        controleDeTensaoImag=cat(2,controleDeTensaoImag,controleDeTensaoRealAux); // Imaginário
        
        // Expansão do vetor para que possa ser inseridas variáveis para cálculo de tensão.
        controleDeTensaoReal=cat(2,controleDeTensaoReal,zeros(size(controleDeTensaoReal,"r"),size(controleDeTensaoReal,"r")+size(controleDeTensaoImag,"r"))); // Real.
        controleDeTensaoImag=cat(2,controleDeTensaoImag,zeros(size(controleDeTensaoImag,"r"),size(controleDeTensaoReal,"r")+size(controleDeTensaoImag,"r"))); // Imagiário.
        
        // Criação do vetor para restição de tensão mínima na barra de geração.
        limiteInferiorParaTensao=zeros(size(controleDeTensaoReal,"r"),size(controleDeTensaoReal,"c"));
        
        // Inserção de variáveis para calculo de tensão.
        controleDeTensaoReal(1,2*(NR+M(1,3))+1)=1; // Real.
        controleDeTensaoImag(1,2*(NR+M(1,3))+2)=1; // Imaginário.
        
        // Identificação de qual variável deverá controlar a tensão.
        limiteInferiorParaTensao(1,2*(NR+M(1,3))+1)=-1;
//    elseif(fonteAMais==2)
//        disp("Você não optou por inserir uma nova fonte de geração");
//        controleDeTensaoReal=[];
//        controleDeTensaoImag=[];
//        limiteInferiorParaTensao=[];
//        M=M(1:NR,:);
//        continue;
//    else
//        disp("Você digitou um valor inválido.");
//    end
endfunction

// ---------------------------------------------------------------
// Funcao responsavel pela construção da Matriz incidência C.
// Entrada : Estrutura do arquivo completo e Matriz laços externos.
// Saída : Matriz incidência C, informação de colunas para resolução
// de caso real, informação de colunas para resolução de caso completo
// de rede.
//----------------------------------------------------------------

function [C,qnt_coluna_MC_incidencia,qnt_coluna_MC,me] = MatrizC(M,controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,fonteAMais)
    
//    if fonteAMais==2 then
//        NR=NR-1;
//        NB=NB-1;
//    end
    AuxC=1;// Auxiliar para criação da matriz inciência de "carga", endereço da coluna.
    AuxG=1;// Auxiliar para criação da matriz inciência de "geração", endereço da coluna.
        
    MAg=zeros(M(1,2),M(1,3)); // Criação de variáveis novas para descartar restrições de desigualdade.

    for i=2:(NR+1)
        // Orientação.
        origem=M(i,1);
        destino=M(i,2);
    
        // Parte real matriz incidencia (A).
        C(origem,AuxC)=1;
        C(destino,AuxC)=-1;
        
        AuxC=AuxC+1;
        
        // Complemento da matriz incidencia (A) com as barras de geração.
        if(M(i,7)==1)
            MAg(origem,AuxG)=1;
            AuxG=AuxG+1;
        end
    end
    
    // Concatenação da Matriz indicencia com variáveis de folga.
    C=cat(2,C,MAg);
    
    // Criação da matriz que será responsável por receber quais ramos devem estar normalmente abertas ou fechadas.
    matrizRestricao=zeros(NR,size(C,"c"));
    
    // Estrutura de repetição responsável pela criação da matriz responsável pela identificação das linhas abertas, assim como ligar ou desligar ramos da rede

    for i=2:(NR+1)

        if M(i,8)==1
            matrizRestricao(i-1,i-1)=1;
        end
    end
    // Ajuste para inserção da variável de folga
    matrizRestricao=cat(1,matrizRestricao,zeros(M(1,3),size(matrizRestricao,"c")));

    // Informação do tamanho da matriz incidência A. Definição para a quantidade de valores presentes na matriz Q.
    qnt_coluna_MC_incidencia=size(C,'c');

    // Criação de matriz auxiliar de mesmo tamanho que a matriz incidencia após a concatenação com variáveis de folga.
    C1=zeros(size(C,'r'),size(C,'c'));
    
    // Criação da parte imaginária da matriz incidencia.
    // A concatenção é feita desta forma para que sejam inseridas novas variáveis, sendo responsáveis pela parte imaginária.
    C2=cat(2,C1,C);
    
    // Junção da matriz incidecia real com seu complemento de zeros.
    // Esta concatenação é feita para que possam ser feita posteriormente a junção com a parte imaginária da rede.
    C=cat(2,C,C1);
    
    // Junção da matriz incidecia real com a parte imaginária.
    C=cat(1,C,C2);

    // Junção de restrição para desligar ramos para partes Real e Imaginária.
    matrizRestricaoAux=matrizRestricao;
    matrizRestricao=cat(2,matrizRestricao,zeros(size(matrizRestricao,"r"),size(matrizRestricao,"c")));
    matrizRestricaoAux=cat(2,zeros(size(matrizRestricaoAux,"r"),size(matrizRestricaoAux,"c")),matrizRestricaoAux);
    matrizRestricao=cat(1,matrizRestricao,matrizRestricaoAux);

    // Junção da matriz Indicência com a matriz restrição.
    C=cat(1,C,matrizRestricao);
//    pause
//    if fonteAMais==1 then
        C=cat(2,C,zeros(size(C,"r"),size(controleDeTensaoImag,"c")-size(C,"c")));
    
        //Inserção da matriz de laço junto à matriz incidência C.
        C=cat(1,C,controleDeTensaoReal);
        C=cat(1,C,controleDeTensaoImag);
        
        //Valor responsável para a criação da matriz P futuramente.
        qnt_coluna_MC=size(C,'c');
        me=size(C,'r');
        //Controle de tensão.
        C=cat(1,C,limiteInferiorParaTensao);
        
//    else
        //Valor responsável para a criação da matriz P futuramente.
//        qnt_coluna_MC=size(C,'c');
//    end
        C=cat(2,C,zeros(size(C,"r"),size(C,"r")-size(C,"c")));
endfunction

// --------------------------------------------------------------------
// Funcao responsavel pela contrução da Matriz de carga b.
// Entrada: Estrutura do arquivo completo, informação de colunas para 
// resolução de caso real, informação de colunas para resolução de caso
// completo de rede.
// Saída : Matriz de carga b.
//---------------------------------------------------------------------

function [b] = MatrizB(M,controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,fonteAMais)
    
//    if fonteAMais==2 then
//        NR=NR-1;
//        NB=NB-1;
//    end
    // Definição inicial dos vetores de carga real e imaginário.
    b=zeros(NB,1); // Real.
    b1=zeros(NB,1); // Imaginário.

    // Criação da matriz B com todos valores negativos.
    for i=2:NR+1
        // Verificação da característica da barra, 1=geração, 
        if(M(i,7)==1)
            // Matriz b parte real da carga.
            b(M(i,1),1)=(M(i,3))//+M((i+1),9));
    
            //Matriz B parte imaginário da carga.
            b1(M(i,1),1)=M(i,4);
        elseif(M(i,3)~=0)
            b(M(i,2),1)=(-1)*(M(i,3));
            
            b1(M(i,2),1)=(-1)*(M(i,4));
        end
    end

    // Complemento referente a barras inseridas de folga e restrições de abertura/fechamento de ramos, além das equações de laço.
//    if(NB==NR) // Fora necessário inserir esta comparação pois no caso da rede de 400barras, haviam na verdade 402 barras e 402 ramos, mas como 1 ramo era feeder e é considerado que neste caso seja inserido uma barra nova, a quantidade de barras e ramos era diferente, ou seja 403 barras e 402 ramos no total.
        // Junção da matriz Solicitação de carga partes Real e Imaginária.
        b=cat(1,b,b1);
        // Inserção de restrições de laço pela lei de Kirchoff para a parte Real da rede.
        b=cat(1,b,zeros(2*(NR+M(1,3)),1));
//        if fonteAMais==1 then
            // Inserção de tensão inicial Real de barra.
            b=cat(1,b,vrReal*ones(size(controleDeTensaoReal,"r"),1));
            
            // Inserção de tensão inicial Imaginária de barra.
            b=cat(1,b,vrImag*ones(size(controleDeTensaoImag,"r"),1));
            
            // Definição do limite inferior para controle de tensão.
            b=cat(1,b,vi*(-1)*ones(size(limiteInferiorParaTensao,"r"),1));
//        end
//    elseif(NB<NR)
//        // Junção da matriz Solicitação de carga partes Real e Imaginária.
//        b=cat(1,b,b1);
//        // Inserção de restrições de laço pela lei de Kirchoff para a parte Real da rede.
//        b=cat(1,b,zeros(2*(NR+M(1,3)),1));
//        
////        if fonteAMais==1 then
//            // Inserção de tensão inicial Real de barra.
//            b=cat(1,b,vrReal*ones(size(controleDeTensaoReal,"r"),1));
//            
//            // Inserção de tensão inicial Imaginária de barra.
//            b=cat(1,b,vrImag*ones(size(controleDeTensaoImag,"r"),1));
//            
//            // Definição do limite inferior para controle de tensão.
//            b=cat(1,b,vi*(-1)*ones(size(limiteInferiorParaTensao,"r"),1));
////            b=cat(1,b,zeros(22,1));
////        end
//    else
//        // Ajuste caso o número de barras for superior ao número de ramos para parte imaginária
////        pause
////        b=cat(1,b,zeros(M(1,3),1));
//        // Inserção das solicitações de carga para parte imaginária da rede.
//        b=cat(1,b,b1);
////        pause
//        // Ajuste caso o número de barras for superior ao número de ramos para parte imaginária
////        b=cat(1,b,zeros(M(1,3),1));
////        pause
//        // Inserção de restrições de laço pela lei de Kirchoff para a parte Real da rede.
//        b=cat(1,b,zeros(2*(NR+M(1,3)),1));        
////        pause
////        if fonteAMais==1 then
//            // Inserção de tensão inicial Real de barra.
//            b=cat(1,b,vrReal*ones(size(controleDeTensaoReal,"r"),1));
//            
//            // Inserção de tensão inicial Imaginária de barra.
//            b=cat(1,b,vrImag*ones(size(controleDeTensaoImag,"r"),1));
//            
//            // Definição do limite inferior para controle de tensão.
//            b=cat(1,b,vi*(-1)*ones(size(limiteInferiorParaTensao,"r"),1));
////        end
//    end
endfunction

//-------------------------------------------------------------------
// Função responsável pela criação da matriz Q.
// Entrada : Estrutura do arquivo completo, informação de colunas para
// resolução de caso real, informação de colunas para resolução de caso
// completo de rede.
// Saída : Matriz simétrica para resistência e reatância Q.
//-------------------------------------------------------------------

function [Q]=MatrizQ(M,qnt_coluna_MC,controleDeTensaoReal,controleDeTensaoImag,C)
    //Dimensionamento de tamanho para matriz Q e auxiliares.
    Q=zeros(NR,qnt_coluna_MC_incidencia);
    Q2=Q;
    Q3=Q;
    
    // Definição das posições onde se encontram barras com 
    numeroDeBarrasComCarga=find(M(2:size(M,"c")-M(1,3),3)~=0);
    
    // Criação da matriz simétrica para resistência e reatância de cada ramo fechado da rede.
    for i=1:(size(numeroDeBarrasComCarga,"c"))
        // Matriz simétrica para resistência.
        Q(i,i)=M(i+1,5);
        
        // Matriz simétrica para reatância.
        Q2(i,i)=M(i+1,6);
    end

    // Criação da matriz simétrica para resistência e reatância de cada ramo aberto da rede.
    for i=(size(numeroDeBarrasComCarga,"c"))+1:NR
        // Matriz simétrica para resistência.
        Q(i,i)=M(i+1,5);
        
        // Matriz simétrica para reatância.
        Q2(i,i)=M(i+1,6);
    end

    // Criação de complemento que será responsável pela resistência e impedância das variáveis de folga.
    AuxQ=zeros(M(1,3),M(1,3));
    AuxQ=cat(2,zeros(M(1,3),NR),AuxQ);

    // Inserção de peso para variáveis de folga.
    Q=cat(1,Q,AuxQ); // Real.

    Q2=cat(1,Q2,AuxQ); // Imaginário.

    //Complemento à matriz Q (dobrando seu tamanho) para que seja possível, posteriormente, fazer a concatenação com a parte referente a reatancias da rede.
    Q=cat(2,Q,zeros(size(Q,'r'),size(Q,'c'))); // Parte real (resistencias).
    Q2=cat(2,zeros(size(Q2,'r'),size(Q2,'c')),Q2); // Parte imaginária (impedâncias).

    // Junção da matriz responsáveis pela resistencia dos ramos da rede e dos pesos para religamento da rede com a matriz responsável pela reatância dos ramos da rede.
    Q=cat(1,Q,Q2);

    // Cálculo de quantos pesos deverão ser adicionados à matriz Q referente às equações de laço.
    quantasBarrasForamAdicionadas=qnt_coluna_MC-size(Q,"c");
    // Expansão para acrescimo de peso para variáveis de laço.
    Q=cat(2,Q,zeros(size(Q,"r"),quantasBarrasForamAdicionadas));
    Q=cat(1,Q,zeros(quantasBarrasForamAdicionadas,size(Q,"c")));

    // Expanção de colunas decorrente da inserção das equações de laço.
    Q=cat(2,Q,zeros(size(Q,"r"),size(C,"r")-size(Q,"c")));
    // Expanção de linhas decorrente da inserção das equações de laço.
    Q=cat(1,Q,zeros(size(C,"r")-size(Q,"r"),size(Q,"c")));

// Correção para que a matriz se torne simétrica, ou seja, diferente de zero diagonal princial.
for i=1:size(Q,"c")
    if Q(i,i)==0
        Q(i,i)=0.0000001;
    end
end

endfunction

//----------------------------------------------------
// Criação matriz p.
// Entrada : Valor de colunas da matriz incidência A.
// Saída : Matriz peso p.
//----------------------------------------------------

function [p]=MatrizP(C)
    
    //Dimensão desta matriz deve ser igual a quantidade de colunas da matriz incidência A.
    a=0;
    for i=1:size(C,"c")
        p(i,1)=0;
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
        disp("Houve um Defeito Falha em algum ramo?");
        algumRamoFalhou=input("Digita 1 (um) para sim ou 0 (zero) para não. ");
        // Verificação se houve Defeito Falha.
        if algumRamoFalhou==0
            break;
        elseif algumRamoFalhou==1 // Confirmação de Defeito Falha.
            //  Comunicação ao usuário.
            defeitoFalha=input("Digite o ramo com Defeito Falha: ");
            // Verificação de possibilidade de existir o Defeito Falha informado.
            if(defeitoFalha>NB)
                disp("Você digitou um valor inválido");
            elseif(restricao~=0) // Confirmação de Defeito Falha informado.
                // Inserção de valor 1 para iniciar uma restrição capaz de desligar o ramo informado tanto para parte real quanto imaginária.
                C(2*NB+defeitoFalha,defeitoFalha)=1;
                C(2*NB+NR+M(1,3)+defeitoFalha,NR+M(1,3)+defeitoFalha)=1;
            elseif(defeitoFalha==0)
                continue;
            end
            
            // Comunicação com o Usuário.
            on_off=input("Deseja ligar alguma linha? 1 para sim, 0 para não. ");
            
            // Confirmação do usuário para ligar algum ramo manualmente.
            if(on_off==1)
                // Demonstração de opções para religamento.
                disp("As opções são:");
                for i=1:NR
                    if(M(i+1,8)==1)
                        a=M(i+1,1);
                        b=M(i+1,2);
                        c="-";
                        disp("Origem Destino - Número da linha");
                        disp(a, b, c, i);
                    end
                end
                
                // Possibiliade de fechamento das opções já fornecidas ao usuário.
                while 1>0
                    // Comunicação ao usuário.
                    ligar=input("Quais linhas que deseja ligar? Quando já colocou todas suas opções, digite 0. ");
                    if(ligar==0) // Foram fechadas todas que o usuário quis.
                        break;
                    elseif(ligar>NR) // Usuário digitou informação inválida.
                        disp("Você digitou um valor inválido")
                        
                    else // Usuário digitou informação válida
                        // Inserção de valor 0 para iniciar uma restrição capaz de ativar o ramo informado tanto para parte real quanto imaginária.
                        C(2*NB+ligar,ligar)=0; // Real.
                        C(2*NB+M(1,3)+NR+ligar,NR+M(1,3)+ligar)=0; // Imaginário.
                    end
                end
                break;
            elseif(on_off==0)
                break;
            else
                disp("Você não digitou um valor válido.");
            end
        else
            disp("Você digitou uma informação inválida");
        end
    end
endfunction

// ------------------------------------------------
// Estutura principal do Algoritimo para Otimização.
// do fluxo de Corrente Alternada com Contingência.
// ------------------------------------------------

// Instrução para criação da matriz responsável pela criação dos laços externos da rede.
[controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,M,fonteAMais]=controleDeTensao(M);

// Instrução para criação da matriz incidência A.
[C,qnt_coluna_MC_incidencia,qnt_coluna_MC,me]=MatrizC(M,controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,fonteAMais);

// Instrução para criação da matriz incidência Q.
[Q]=MatrizQ(M,qnt_coluna_MC,controleDeTensaoReal,controleDeTensaoImag,C);

// Instrução para criação da matriz incidência P.
[p]=MatrizP(C);

//[C]=restricao(C,M);

// Instrução para criação da matriz incidência B.
[b]=MatrizB(M,controleDeTensaoReal,controleDeTensaoImag,limiteInferiorParaTensao,fonteAMais);

// Limite inferior.
ci=(-5)*ones(size(C,"c"),1);
// Limite superior.
cs=ci*(-1);

// Função de otimização QPSOLVE - objetivo: Minimização de perdas na rede.
[xopt,iact,iter,fopt]=qpsolve(Q,p,C,b,ci,cs,me)
//
// Informação ao usuário da resolução da rede.
if xopt~=[]
    disp("Solução Ótima Encontrada!");
    disp(fopt,"O valor ótimo encontrado para a função objetivo.");
    disp(iter,"Foram necessárias esta quantidade de iterações para convergir.");
    disp("O primeiro valor se refere às iterações e o segundo às restrições desativadas para resolução da rede.");
else
    disp("Solução não encontrada.");
end

// Variáveis de resolução de rede.
xopt=xopt(1:((2*(NR+M(1,3))+(size(controleDeTensaoReal,'r')+size(controleDeTensaoImag,'r'))),1));
