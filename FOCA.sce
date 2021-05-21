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
entradaDeDados=input("Digite o endereço: ");//com localização do arquivo de entrada de dados. Obs: seguir instruções no arquivo Guideline.

// Arquivo de Entrada com a estrutura da Rede de Distribuição.
M = fscanfMat(entradaDeDados, "%lg"); // Importa arquivo que apresenta os dados de entrada da rede. 

//"C:\Users\Higor\Desktop\Iniciação científica\2021\Rede 5 barras\rede5barras.txt"
//"C:\Users\Higor\Desktop\Iniciação científica\2021\Rede 16 barras\rede16barras.txt"
//"C:\Users\Higor\Desktop\Iniciação científica\2021\Rede Nova\reduzida\rede28barras.txt" ------ rede reduzida

//"C:\Users\Higor\Desktop\Iniciação científica\2021\Rede 34 barras\rede34barras.txt" -------- rede 34 barras dissertação Adolfo
//"C:\Users\Higor\Desktop\Iniciação científica\2021\Rede 34 barras\rede34barrasAdmitancia.txt"
//"C:\Users\Higor\Desktop\Iniciação científica\2021\400barras\rede400barras.txt"
//"C:\Users\Higor\Desktop\Iniciação científica\2021\400barras\400barrasV2.txt" -------- rede 400 barras dissertação Adolfo

// Salva o número de ramos
NR=M(1,1);

// Salva o número de barras
NB=M(1,2)+M(1,3); 

//-----------------------
// Criação matriz c
// Matriz responsável pela criação dos caminhos da rede, relacionado aos seus ramos.
//-----------------------
function [C,qnt_coluna_MA_incidencia,qnt_coluna_MA] = MatrizA(M,LE,quantosLacos)

    AuxC=1;// Auxiliar para criação da matriz inciência de "carga", endereço da coluna.
    AuxG=1;// Auxiliar para criação da matriz inciência de "geração", endereço da coluna.
        
    MAg=zeros(M(1,2),M(1,3)); //Criação de variáveis novas para descartar restrições de desigualdade.
        
    for i=2:(NR+1)
        //Orientação.
        origem=M(i,1);
        destino=M(i,2);
    
        //Parte real matriz incidencia (A).
        C(origem,AuxC)=1;
        C(destino,AuxC)=-1;
        
        AuxC=AuxC+1;
        
        //Complemento da matriz incidencia (A) com as barras de geração.
        if(M(i,7)==1)
            MAg(origem,AuxG)=1;
            AuxG=AuxG+1;
        end
    end
    
    //Concatenação da Matriz indicencia com variáveis de folga.
    C=cat(2,C,MAg);
    
    //Criação da matriz que será responsável por receber quais ramos devem estar normalmente abertas ou fechadas.
    matrizRestricao=zeros(NB,size(C,"c"));
    
    // Estrutura de repetição responsável pela criação da matriz responsável pela identificação das linhas abertas, assim como ligar ou desligar ramos da rede
    for i=2:(NR+1)
        if M(i,8)==1
            matrizRestricao(i-1,i-1)=1;
        end
    end
    //Ajuste para inserção da variável de folga
    matrizRestricao=cat(1,matrizRestricao,zeros(M(1,3),size(matrizRestricao,"c")));

    // Informação do tamanho da matriz incidência A. Definição para a quantidade de valores presentes na matriz Q.
    qnt_coluna_MA_incidencia=size(C,'c');
    
    // Criação de matriz auxiliar de mesmo tamanho que a matriz incidencia após a concatenação com variáveis de folga.
    C1=zeros(size(C,'r'),size(C,'c'));
    
    // Criação da parte imaginária da matriz incidencia.
    // A concatenção é feita desta forma para que sejam inseridas novas variáveis, sendo responsáveis pela parte imaginária.
    C2=cat(2,C1,C);
    
    // Junção da matriz incidecia real com seu complemento de zeros.
    // Esta concatenação é feita para que possam ser feita posteriormente a junção com a parte imaginária da rede.
    C=cat(2,C,C1);
    
    //Junção da matriz incidecia real com a parte imaginária.
    C=cat(1,C,C2);
    
    //Ajuste para desligamento de ramo
    C=cat(2,C,zeros(size(C,"r"),size(C,"c")));

    //Ajuste para concatenação das equações de laço junto à matriz incidência C.
    LE=cat(2,LE,zeros(size(LE,'r'),size(C,"c")-size(LE,'c')));

    // Comando para junção vertical entre matrizes
    matrizRestricaoAux=matrizRestricao;
    matrizRestricao=cat(2,matrizRestricao,zeros(size(matrizRestricao,"r"),size(matrizRestricao,"c")));
    matrizRestricaoAux=cat(2,zeros(size(matrizRestricaoAux,"r"),size(matrizRestricaoAux,"c")),matrizRestricaoAux);
    matrizRestricao=cat(1,matrizRestricao,matrizRestricaoAux);
    matrizRestricao=cat(2,matrizRestricao,zeros(size(matrizRestricao,"r"),size(matrizRestricao,"c")));

    C=cat(1,C,matrizRestricao);

    //Inserção da matriz de laço junto à matriz incidência C.
    C=cat(1,C,LE);

    C=cat(2,C,zeros(size(C,"r"),quantosLacos));
    
    //Valor responsável para a criação da matriz P futuramente.
    qnt_coluna_MA=size(C,'c');
endfunction

//-----------------------
// Criação matriz b.
// Esta matriz é responsável pela definição das cargas e gerações da rede.
//-----------------------

function [b] = MatrizB(M,quantosLacos,qnt_coluna_MA_incidencia)
    // Criação da matriz B com todos valores negativos.
    for i=1:NR
        // Matriz b parte real da carga.
        b(i)=(-1)*(M(i+1,3))//+M((i+1),9));

        //Matriz B parte imaginário da carga.
        b1(i)=(-1)*M(i+1,4);
    end
    //Ajuste para casos onde o número de barras é maior que o número de ramos.
    if NB>NR then
        b=cat(1,b,zeros(M(1,3),1));
        b1=cat(1,b1,zeros(M(1,3),1));
    end
    // Atualização de sinais da matriz B.
    // Sendo positivo o que é gerado e negativo o que é consumido.
    for i=1:NR
        if(M(i+1,7)==1)
            //Correção para a parte Real.
            origem=M(i+1,1);
            b(origem)=b(origem)*(-1);
            
            //Correção para a parte imaginária.
            b1(origem)=b1(origem)*(-1);
        end
    end
    
    //Complemento referente a barras inseridas de 'sobra' da rede tanto em real, quanto imaginário.
    if(NB==NR)// Fora necessário inserir esta comparação pois no caso da rede de 400barras, haviam na verdade 402 barras e 402 ramos, mas como 1 ramo era feeder e é considerado que neste caso seja inserido uma barra nova, a quantidade de barras e ramos era diferente, ou seja 403 barras e 402 ramos no total.
        // Ajuste de tamanho para matriz referente a Solicitação de Carga parte Real
//        b=b(1:NB)
//        // Ajuste caso o número de barras for superior ao número de ramos para parte Real
//        b=cat(1,b,zeros(M(1,3),1));
//        // Ajuste de tamanho para matriz referente a Solicitação de Carga parte Imaginária
//        b1=b1(1:NB);
//        // Ajuste caso o número de barras for superior ao número de ramos para parte Real
//        b1=cat(1,b1,zeros(M(1,3),1));
        // Junção da matriz Solicitação de carga partes Real e Imaginária.
        b=cat(1,b,b1);
//        pause
        // Inserção de restrições de laço pela lei de Kirchoff para a parte Real da rede.
        b=cat(1,b,zeros(2*qnt_coluna_MA_incidencia,1));
//        b=cat(1,b,zeros(quantosLacos,1));
//pause

        // Inserção de restrições de laço + desligamento de ramos pela lei de Kirchoff para a parte Imaginária da rede.
//        b=cat(1,b,zeros(size(b1,"r")+quantosLacos,1));
//        pause
    b=cat(1,b,ones(quantosLacos,1));
//    C:\Users\Higor\Documents\OtimizacaoDeRedes\ArquivosDeEntradaDeDados\FOCA
//pause
    elseif(NB<NR)
        // Ajuste de tamanho para matriz referente a Solicitação de Carga parte Real e Imaginária
        b=b(1:NB);// Real
        b1=b1(1:NB);// Imaginária
        // Junção da matriz Solicitação de carga partes Real e Imaginária.
        b=cat(1,b,b1);
        // Inserção de restrições para abrir ou fechar um determinado arco/ramo
        b=cat(1,b,zeros(2*qnt_coluna_MA_incidencia,1));//b=cat(1,b,zeros(NR,1));
        // Inserção de restrições de laço  pela lei de Kirchoff para a parte Imaginária da rede.
        b=cat(1,b,ones(quantosLacos,1));//b=cat(1,b,zeros(M(1,3)+quantosLacos,1));
        
    else
        // Ajuste de tamanho para matriz referente a Solicitação de Carga parte Real e Imaginária
        b=b(1:NB);
        b1=b1(1:NB);

        //Ajuste caso o número de barras for superior ao número de ramos para parte imaginária
        b=cat(1,b,zeros(M(1,3),1));
        // Inserção das solicitações de carga para parte imaginária da rede.
        b=cat(1,b,b1);
        
        //Ajuste caso o número de barras for superior ao número de ramos para parte imaginária
        
        b=cat(1,b,zeros(M(1,3),1));
        
        //Ajuste caso o número de barras for superior ao número de ramos para parte real
        
        // Inserção de restrições de laço pela lei de Kirchoff para a parte Real da rede.
        b=cat(1,b,zeros(2*qnt_coluna_MA_incidencia,1));//b=cat(1,b,zeros(NR,1));
        
        // Inserção de restrições de laço  pela lei de Kirchoff para a parte Imaginária da rede.
        b=cat(1,b,ones(quantosLacos,1));//b=cat(1,b,zeros(M(1,3)+quantosLacos,1));
    end
endfunction

//-----------------------
//Criação matriz Q.
//-----------------------

function [Q]=MatrizQ(M,qnt_coluna_MA_incidencia,quantosLacos)
    
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

    //Incremento da matriz auxiliar Q1 à matriz principal real e imaginária.
    Q=cat(1,Q,Q1); //Real.
    Q2=cat(1,Q2,Q1); //Imaginária.
//    pause
    //Complemento à matriz Q (dobrando seu tamanho) para que seja possível, posteriormente, fazer a concatenação com a parte referente a reatancias da rede
    Q=cat(2,Q,zeros(size(Q,'r'),size(Q,'c')));

    //------------------------------------------
    //teste para incrementar desligamento de ramos
    Qx=zeros(qnt_coluna_MA_incidencia,qnt_coluna_MA_incidencia); //Organização 
    Qy=0.0000001*eye(qnt_coluna_MA_incidencia,qnt_coluna_MA_incidencia);
    Q3=cat(2,Qx,Qy);
    Q=cat(1,Q,Q3);
    Q=cat(2,Q,zeros(size(Q,'r'),size(Q,"c")));
    
    Qx=cat(2,Qx,Qx);
    Q3=cat(2,Qx,Q3);
    
    Q2=cat(2,zeros(qnt_coluna_MA_incidencia,qnt_coluna_MA_incidencia),Q2);
    Q2=cat(2,zeros(qnt_coluna_MA_incidencia,qnt_coluna_MA_incidencia),Q2);
    Q2=cat(2,Q2,zeros(qnt_coluna_MA_incidencia,qnt_coluna_MA_incidencia));
    
    Q=cat(1,Q,Q2);
    Q=cat(1,Q,Q3);
    
    
    //teste
    Q=cat(2,Q,zeros(size(Q,"r"),quantosLacos));
    Q7=zeros(quantosLacos,size(Q,"c")-quantosLacos);
    Q7=cat(2,Q7,0.0000001*eye(quantosLacos,quantosLacos));
    Q=cat(1,Q,Q7);

endfunction

//-----------------------
//Criação matriz p
//-----------------------
function [p]=MatrizP(qnt_coluna_MA)
    
    //Dimensão desta matriz deve ser igual a quantidade de colunas da matriz incidência A.
    a=0;
    for i=1:qnt_coluna_MA
        p(i,1)=0;
    end

endfunction

function [C]=restricao(C,M)
    while 1>0 do
        restricao=input("Digite o ramo com Defeito Falha: ");
        if(restricao>NB)
            disp("Você digitou um valor inválido");
        elseif(restricao~=0)
            C(2*NB+restricao,restricao)=1;
            C(2*NB+M(1,3)+NR+restricao,NR+restricao+M(1,3))=1;
            disp(2*NB+M(1,3)+NR+restricao,NR+restricao+M(1,3));
        elseif(restricao==0)
            continue;
        end
        
        on_off=input("Deseja ligar alguma linha? 1 para sim, 0 para não. ")
        if(on_off==1)
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
            while 1>0
                ligar=input("Quais linhas que deseja ligar? Quando já colocou todas suas opções, digite 0. ");
                if(ligar==0)
                    break;
                elseif(ligar>NR)
                    disp("Você digitou um valor inválido")
                else
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
    end
endfunction

function [LE,quantosLacos]=lacosExternos(M)
    mAux=M(2:NR+1,1:size(M,'c'));
    // -----------------------------------------------------
    // Funcao responsável pela construçao do laço externo 
    // Entrada : Estrutura do arquivo completo 
    // Saida : Matriz de lacos externos
    // -----------------------------------------------------
    
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
    
            //Busca amontante até feeder, ou seja, barra de geração
            while (1>0)
                origem=find(D(:,1)==buss); //Mostra a posição onde se encontra o próximo ramo para formação de laços.
                laco(origem(1,1),1)=mAux(origem(1,1),5); //Monta o vetor caminhos, mostrando os laços feitos.
                buss=O(origem(1,1),1); //Recebe o valor que se encontra na posição acima destacada dentro do vetor origem.
                parada=find(O(:,1)==buss(1,1)) //Procura onde está a origem do buss, retornando a posição da matriz O.
                parou=mAux(parada,1) //Atribui o valor referente a posição da matriz O encontrada na linha superior.
                condicaoDeParada=find(geradores(:,:)==parou) //procura na matriz geradores, se já chegou em algum deles
                
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
    
    quaisLacos2=quaisLacos;
    
    barraLaco=cat(1,barraLaco,barraLacoImag);
    
    quaisLacos=cat(2,quaisLacos,quaisLacosImag);
    
    quaisLacosImag=cat(2,quaisLacosImag*(-1),quaisLacos2);
    
    quaisLacos=cat(1,quaisLacos,quaisLacosImag);
    
    quaisLacos=quaisLacos'
    
    quaisLacos=cat(2,quaisLacos,eye(size(quaisLacos,'r'),size(quaisLacos,'r')));
    LE=quaisLacos;
    quantosLacos=size(LE,'r');
endfunction

//Instrução para criação da matriz responsável pela criação dos laços externos da rede.
[LE,quantosLacos]=lacosExternos(M);

//Instrução para criação da matriz incidência A.
[C,qnt_coluna_MA_incidencia,qnt_coluna_MA]=MatrizA(M,LE,quantosLacos);

//Instrução para criação da matriz incidência Q.
[Q]=MatrizQ(M,qnt_coluna_MA_incidencia,quantosLacos);

//Instrução para criação da matriz incidência P.
[p]=MatrizP(qnt_coluna_MA);

//[C]=restricao(C,M);

//Instrução para criação da matriz incidência B.
[b]=MatrizB(M,quantosLacos,qnt_coluna_MA_incidencia);

ci=[];
cs=[];
//ci=(-5)*ones(size(C,"c"),1);
//cs=ci*(-1);

[xopt,iact,iter,fopt]=qpsolve(Q,p,C,b,ci,cs,qnt_coluna_MA)

xopt=xopt(1:((2*(NR+M(1,3))+(size(LE,'r'))),1));

if xopt~=[]
    disp("Solução Ótima Encontrada!")
    disp(fopt,"O valor ótimo encontrado para a função objetivo.")
else
    disp("Solução não encontrada.")
end

