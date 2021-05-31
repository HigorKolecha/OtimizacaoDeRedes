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

// Solicitação ao usuário do endereço para obtenção do arquivo de entrada.
entradaDeDados=input("Digite o endereço com localização do arquivo de entrada de dados. Obs: seguir instruções no arquivo Guideline. ");

// Arquivo de Entrada com a estrutura da Rede de Distribuição.
M = fscanfMat(entradaDeDados, "%lg"); // Importa arquivo que apresenta os dados de entrada da rede. 

// Inicio de contagem de timer para tempo de resolução da rede.
tic();

// Salva o número de ramos.
NR=M(1,1);
// Salva o número de barras.
NB=M(1,2)+M(1,3); 
// Valor de refencia de tensão inicial Real.
vf=1.0;
// Valor de refencia de tensão inicial Imaginário.
vfi=0;

// --------------------------------------------------------------------
// Funcao responsável pela construçao do laço externo.
// Entrada : Estrutura do arquivo completo.
// Saida : Matriz de lacos externos, informação de quantos laços foram
// criados.
// -------------------------------------------------------------------

function [LE,quantosLacos,barraLaco,barraLacoImag]=lacosExternos(M)
    // Ajuste da variável M para retirada de informações gerais.
    mAux=M(2:NR+1,1:size(M,'c'));
    
    // Criação de vetores auxiliares
    O=mAux(:,1); // Criação do vetor origem.
    D=mAux(:,2); // Criação do vetor destino.
    
    // Criação do caminho que a corrente faz para poder criar um laço para a parte real e imaginária da rede.
    quaisLacos=[]; // Real.
    quaisLacosImag=[];// Imaginária.
    
    // Quais barras tiveram seus laços criados para a parte real e imaginária da rede. 
    barraLaco=[]; // Real.
    barraLacoImag=[]; // Imaginária.
    
    // Definição de posições das barras geradoras na matriz mAux.
    geradores=find(mAux(:,7)==1);
    
    // Criação uma matriz com todas barras geradoras da rede.
    geradores=mAux(geradores,1);
    
    // Estrutura de repetição para criação de laços da rede.
    for i=NR:-1:1
        // Definição da barra inicial que posteriormente será verificada a existência de laços.
        buss=mAux(i,2);
        
        // Laços para parte real da rede.
        if (mAux(i,3)~=0) // Verificação da existencia de carga.
            laco=zeros(NR+M(1,3),1); // Criação da matriz referente ao caminho que a corrente faz para poder criar um laço.
            
            barraLaco=cat(1,barraLaco,i); // Faz a concatenação da matriz "barralaco" a fim de serem adicionadas apenas as barras com laços existentes.
            
            // Busca amontante até feeder.
            while (1>0)
                origem=find(D(:,1)==buss); // Mostra as possíveis posições onde se encontra o próximo ramo para formação de laços.
                
                laco(origem(1,1),1)=mAux(origem(1,1),5); // Monta o vetor caminhos, mostrando os laços feitos.
                
                buss=O(origem(1,1),1); // Recebe o valor que se encontra na posição acima destacada dentro do vetor origem.
                
                parada=find(O(:,1)==buss(1,1)); // Procura onde está a origem do buss, retornando a posição da matriz O.
                
                parou=mAux(parada,1); // Atribui o valor referente a posição da matriz O encontrada na linha superior.
                
                condicaoDeParada=find(geradores(:,:)==parou); // Procura na matriz geradores, se já chegou em algum deles.
                
                if(condicaoDeParada~=[]) // Caso o valor da condicaoDeParada seja diferente de vazio, quer dizer que ainda não estamos em um gerador, desta forma, a função while deverá continuar. Caso contrário, deverá parar.
                    break
                end
            end
            buss=mAux(i,2); // Retorna o valor inicial referente a barra a ser utilizada neste momento, para que possa ser criados outros laços futuramente.
            quaisLacos=cat(2,quaisLacos,laco); // Concatenação horizontal, inserindo todos laços 'reais' criados em apenas uma matriz.
        end
        
        // Laços para parte imaginária da rede.
        if (mAux(i,4)~=0) // Verificação da existencia de carga. 
            lacoImag=zeros(NR+M(1,3),1); // Criação do caminho que a corrente faz para poder criar um laço.
            barraLacoImag=cat(1,barraLacoImag,i); // Faz a concatenação da matriz "barralaco" a fim de serem adicionadas apenas as barras com laços existentes.

            // Busca amontante até feeder, ou seja, barra de geração.
            while (1>0)
                origem=find(D(:,1)==buss); // Mostra a posição onde se encontra o próximo ramo para formação de laços.
                
                lacoImag(origem(1,1),1)=mAux(origem(1,1),6); // Monta o vetor caminhos, mostrando os laços feitos.
                
                buss=O(origem(1,1),1); // Recebe o valor que se encontra na posição acima destacada dentro do vetor origem.
                
                parada=find(O(:,1)==buss(1,1)) // Procura onde está a origem do buss, retornando a posição da matriz O.
                
                parou=mAux(parada,1); // Atribui o valor referente a posição da matriz O encontrada na linha superior.
                
                condicaoDeParada=find(geradores(:,:)==parou); // Procura na matriz geradores, se já chegou em algum deles.
                if(condicaoDeParada~=[]) // Caso o valor da condicaoDeParada seja diferente de vazio, quer dizer que ainda não estamos em um gerador, desta forma, a função while deverá continuar. Caso contrário, deverá parar.
                    break
                end
            end
            buss=mAux(i,2); // Retorna o valor inicial referente a barra a ser utilizada neste momento, para que possa ser criados outros laços futuramente.
            quaisLacosImag=cat(2,quaisLacosImag,lacoImag); // Concatenação horizontal, inserindo todos laços 'imaginários' criados em apenas uma matriz.
        end
    end
    
    // Criação da matriz auxiliar.
    quaisLacos2=quaisLacos;
    
    // Junção das informações de quais barras geraram laço para parte real e imaginária.
//    barraLaco=cat(1,barraLaco,barraLacoImag);

    // Junção de laços reais e imaginários.
    quaisLacos=cat(2,quaisLacos,quaisLacosImag);
    
    // Junção dos laços para calculo de tensão imaginária.
    quaisLacosImag=cat(2,quaisLacosImag*(-1),quaisLacos2);
    
    // Junção de todos laços reais e imaginários para calculo da tensão de barra parte real e imaginária.
    quaisLacos=cat(1,quaisLacos,quaisLacosImag);
    
    // Transposição da matriz para concatenação com a matriz indicencia posteriormente.
    quaisLacos=quaisLacos';
    
    // Expansão da matriz para resolução da rede.
    quaisLacos=cat(2,quaisLacos,eye(size(quaisLacos,'r'),size(quaisLacos,'r')));
    
    // Nomeação de nova variável, recebendo a matriz de laços criada.
    LE=quaisLacos;
    
    // Dimensão de quantos laços ao todo foram criados.
    quantosLacos=size(LE,'r');
endfunction

// ---------------------------------------------------------------
// Funcao responsavel pela construção da Matriz incidência C.
// Entrada : Estrutura do arquivo completo e Matriz laços externos.
// Saída : Matriz incidência C, informação de colunas para resolução
// de caso real, informação de colunas para resolução de caso completo
// de rede.
//----------------------------------------------------------------

function [C,qnt_coluna_MC_incidencia,qnt_coluna_MC] = MatrizC(M,LE,quantosLacos)
    AuxC=1; // Auxiliar para criação da matriz inciência de "carga", endereço da coluna.
    AuxG=1; // Auxiliar para criação da matriz inciência de "geração", endereço da coluna.
        
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
    matrizRestricao=zeros(NB,size(C,"c"));
    
    // Estrutura de repetição responsável pela criação da matriz responsável pela identificação das linhas abertas, assim como ligar ou desligar ramos da rede.
    for i=2:(NR+1)
        if M(i,8)==1
            matrizRestricao(i-1,i-1)=1;
        end
    end
    
    // Ajuste para inserção da variável de folga.
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
    
    // Ajuste para desligamento de ramo.
    C=cat(2,C,zeros(size(C,"r"),size(LE,"c")-size(C,"c")));

    // Junção de restrição para desligar ramos para partes Real e Imaginária.
    matrizRestricaoAux=matrizRestricao;
    matrizRestricao=cat(2,matrizRestricao,zeros(size(matrizRestricao,"r"),size(matrizRestricao,"c")));
    matrizRestricaoAux=cat(2,zeros(size(matrizRestricaoAux,"r"),size(matrizRestricaoAux,"c")),matrizRestricaoAux);
    matrizRestricao=cat(1,matrizRestricao,matrizRestricaoAux);
//    matrizRestricao=cat(2,matrizRestricao,zeros(size(matrizRestricao,"r"),size(matrizRestricao,"c")));
    
    matrizRestricao=cat(2,matrizRestricao,zeros(size(matrizRestricao,"r"),size(LE,"c")-size(matrizRestricao,"c")));
    
    // Junção da matriz Indicência com a matriz restrição.
    C=cat(1,C,matrizRestricao);

    // Inserção da matriz de laço junto à matriz incidência C.
    C=cat(1,C,LE);
    
    // Inserção de novas colunas para que as restrições de laço funcionem.
    C=cat(2,C,zeros(size(C,"r"),size(C,"r")-size(C,"c")));
    
    //Valor responsável para a criação da matriz P futuramente.
    qnt_coluna_MC=size(C,'c');
endfunction

// --------------------------------------------------------------------
// Funcao responsavel pela contrução da Matriz de carga b.
// Entrada: Estrutura do arquivo completo, informação de colunas para 
// resolução de caso real, informação de colunas para resolução de caso
// completo de rede.
// Saída : Matriz de carga b.
//---------------------------------------------------------------------

function [b] = MatrizB(M,quantosLacos,qnt_coluna_MC_incidencia,barraLaco,barraLacoImag)
    // Criação da matriz B com todos valores negativos.
    for i=1:NR
        // Matriz b parte real da carga.
        b(i)=(-1)*(M(i+1,3))//+M((i+1),9));

        // Matriz B parte imaginário da carga.
        b1(i)=(-1)*M(i+1,4);
    end
    
    // Ajuste para casos onde o número de barras é maior que o número de ramos.
    if NB>NR then
        b=cat(1,b,zeros(M(1,3),1));
        b1=cat(1,b1,zeros(M(1,3),1));
    end
    
    // Atualização de sinais da matriz B.
    // Sendo positivo o que é gerado e negativo o que é consumido.
    for i=1:NR
        if(M(i+1,7)==1)
            // Correção para a parte Real.
            origem=M(i+1,1);
            b(origem)=b(origem)*(-1);
            
            // Correção para a parte imaginária.
            b1(origem)=b1(origem)*(-1);
        end
    end
    
    // Complemento referente a barras inseridas de folga e restrições de abertura/fechamento de ramos.
    if(NB==NR) // Fora necessário inserir esta comparação pois no caso da rede de 400barras, haviam na verdade 402 barras e 402 ramos, mas como 1 ramo era feeder e é considerado que neste caso seja inserido uma barra nova, a quantidade de barras e ramos era diferente, ou seja 403 barras e 402 ramos no total.
        // Junção da matriz Solicitação de carga partes Real e Imaginária.
        b=cat(1,b,b1);
        
        // Inserção de restrições de laço pela lei de Kirchoff para a parte Real da rede.
        b=cat(1,b,zeros(2*(NR+M(1,3)),1));
        
        // Inserção de tensão inicial de barra.
        b=cat(1,b,vf*ones(size(barraLaco,"r"),1));
        b=cat(1,b,vfi*ones(size(barraLacoImag,"r"),1));
    elseif(NB<NR)
        // Ajuste de tamanho para matriz referente a Solicitação de Carga parte Real e Imaginária
        b=b(1:NB);// Real
        b1=b1(1:NB);// Imaginária
        
        // Junção da matriz Solicitação de carga partes Real e Imaginária.
        b=cat(1,b,b1);
        
        // Inserção de restrições de laço pela lei de Kirchoff para a parte Real da rede.
        b=cat(1,b,zeros(2*(NR+M(1,3)),1));
        
        // Inserção de tensão inicial de barra.
        b=cat(1,b,vf*ones(size(barraLaco,"r"),1));
        b=cat(1,b,vfi*ones(size(barraLacoImag,"r"),1));
//        pause
    else
        // Ajuste de tamanho para matriz referente a Solicitação de Carga parte Real e Imaginária.
        b=b(1:NB);// Real.
        b1=b1(1:NB);// Imaginária.

        // Ajuste caso o número de barras for superior ao número de ramos para parte imaginária
        b=cat(1,b,zeros(M(1,3),1));
        
        // Inserção das solicitações de carga para parte imaginária da rede.
        b=cat(1,b,b1);
        
        // Ajuste caso o número de barras for superior ao número de ramos para parte imaginária
        b=cat(1,b,zeros(M(1,3),1));
        
        // Inserção de restrições de laço pela lei de Kirchoff para a parte Real da rede.
        b=cat(1,b,zeros(2*(NR+M(1,3)),1));
        
        // Inserção de tensão inicial de barra.
        b=cat(1,b,vf*ones(size(barraLaco,"r"),1));
        b=cat(1,b,vfi*ones(size(barraLacoImag,"r"),1));
    end
endfunction

//-------------------------------------------------------------------
// Função responsável pela criação da matriz Q.
// Entrada : Estrutura do arquivo completo, informação de colunas para
// resolução de caso real, informação de colunas para resolução de caso
// completo de rede.
// Saída : Matriz simétrica para resistência e reatância Q.
//-------------------------------------------------------------------

function [Q]=MatrizQ(M,qnt_coluna_MC_incidencia,quantosLacos)
    // Dimensionamento de tamanho para matriz Q.
    Q=zeros(NR,qnt_coluna_MC_incidencia);
    Q2=zeros(NR,qnt_coluna_MC_incidencia);
    Q3=zeros(NR,qnt_coluna_MC_incidencia);
    
    // Criação da matriz simétrica para resistência e reatância de cada barra.
    for i=2:(NR+1)
        // Matriz simétrica para resistência.
        Q(i-1,i-1)=M(i,5);
        
        // Matriz simétrica para reatância.
        Q2(i-1,i-1)=M(i,6);
    end
    
    // Criação de complemento que será responsável pela resistência e impedância das as barras de geração.
    AuxQ=0.0000001*eye(M(1,3),M(1,3));
    
    // Responsável pela criação de uma matriz complementar, irá possibilitar a concatenação posteriormente das variáveis de folga da rede.
    Q1=zeros(M(1,3),NR);
    
    // Concatenação de incremento para as variáveis com baixa resistencia referentes às variáveis de folga.
    Q1=cat(2,Q1,AuxQ);

    // Incremento da matriz auxiliar Q1 à matriz principal real e imaginária.
    Q=cat(1,Q,Q1); //Real.
    Q2=cat(1,Q2,Q1); //Imaginária.
    
    // Complemento à matriz Q (dobrando seu tamanho) para que seja possível, posteriormente, fazer a concatenação com a parte referente a reatancias da rede.
    // Parte real (resistencias)
    Q=cat(2,Q,zeros(size(Q,'r'),size(Q,'c')));
    
    // Complemento à matriz Q (dobrando seu tamanho) para que seja possível, posteriormente, fazer a concatenação com a parte referente a reatancias da rede.
    Q2=cat(2,zeros(size(Q2,"r"),size(Q2,"c")),Q2);
    
    // Junção da matriz responsáveis pela resistencia dos ramos da rede e dos pesos para religamento da rede com a matriz responsável pela reatância dos ramos da rede.
    Q=cat(1,Q,Q2);
    
    // Expanção de colunas decorrente da inserção das equações de laço.
    Q=cat(2,Q,zeros(size(Q,"r"),size(C,"c")-size(Q,"c")));
    
    // Expanção de linhas decorrente da inserção das equações de laço.
    Q=cat(1,Q,zeros(size(C,"c")-size(Q,"r"),size(Q,"c")));

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

function [p]=MatrizP(qnt_coluna_MC)
    //Dimensão desta matriz deve ser igual a quantidade de colunas da matriz incidência A.
    a=0;
    for i=1:qnt_coluna_MC
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
        
        // Verificação se houve Defeito Falha
        if algumRamoFalhou==0
            break;
        elseif algumRamoFalhou==1
            // Comunicação ao usuário.
            defeitoFalha=input("Digite o ramo com Defeito Falha: ");
            
            // Verificação de possibilidade de existir o Defeito Falha informado.
            if(defeitoFalha>NB)
                disp("Você digitou um valor inválido");
            elseif(defeitoFalha~=0)
                // Inserção de valor 1 para iniciar uma restrição capaz de desligar o ramo informado tanto para parte real quanto imaginária.
                C(2*NB+M(1,3)+defeitoFalha,M(1,3)+defeitoFalha)=1;
                C(2*NB+NR+M(1,3)+defeitoFalha,NR+defeitoFalha+M(1,3))=1;
            elseif(defeitoFalha==0)
                continue;
            end
            
            // Comunicação com o usuário.
            on_off=input("Deseja ligar alguma linha? 1 para sim, 0 para não. ")
            
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
                    ligar=input("Quais linhas que deseja ligar? Quando já colocou todas suas opções, digite 0. ");
                    if(ligar==0) // Foram fechadas todas que o usuário quis.
                        break;
                    elseif(ligar>NR) // Usuário digitou informação inválida.
                        disp("Você digitou um valor inválido")
                    else // Usuário digitou informação válida.
                        // Inserção de valor 0 para iniciar uma restrição capaz de ativar o ramo informado tanto para parte real quanto imaginária.
                        C(2*NB+ligar,ligar)=0;
                        C(2*NB+M(1,3)+NR+ligar,NR+ligar+M(1,3))=0;
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
[LE,quantosLacos,barraLaco,barraLacoImag]=lacosExternos(M);

// Instrução para criação da matriz incidência A.
[C,qnt_coluna_MC_incidencia,qnt_coluna_MC]=MatrizC(M,LE,quantosLacos);

// Instrução para criação da matriz incidência Q.
[Q]=MatrizQ(M,qnt_coluna_MC_incidencia,quantosLacos);

// Instrução para criação da matriz incidência P.
[p]=MatrizP(qnt_coluna_MC);

// Instrução para inserção de desligamento de barras em casa de Defeito Falha.
[C]=restricao(C,M);

// Instrução para criação da matriz incidência B.
[b]=MatrizB(M,quantosLacos,qnt_coluna_MC_incidencia,barraLaco,barraLacoImag);

// Limite inferior.
ci=(-5)*ones(size(C,"c"),1);
// Limite superior.
cs=ci*(-1);

// Função de otimização QPSOLVE - objetivo: Minimização de perdas na rede.
[xopt,iact,iter,fopt]=qpsolve(Q,p,C,b,ci,cs,qnt_coluna_MC);

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
    disp("Solução não encontrada.")
end

// Variáveis de resolução de rede.
xopt=xopt(1:((2*(NR+M(1,3))+(size(LE,'r'))),1));
