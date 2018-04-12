#!/bin/bash

# Quantidades de elementos:
nelDir1p=200 # Valor padrao de elementos em x
nelDir1=${1:-${nelDir1p}} # Elementos em x
nelDir2p=200 # Valor padrao de elementos em y
nelDir2=${2:-${nelDir2p}} # Elementos em y
Elements=$(( nelDir1 * nelDir2)) # Numero total de elementos  
domI_X=0.0 # Ponto inferior do dominio em X
domF_X=1.0 # Ponto superior do dominio em X
domI_Y=0.0 # Ponto inferior do dominio em Y
domF_Y=1.0 # Ponto superior do dominio em Y

# Instrucao 1:
iecho=0
printf "%5s \n" ${iecho} 

# Instrucao 2:
printf "Poisson Problem - Convergencia: Q2Q2-2 - Refinamento h: teta=-1 (otimo) \n"

# Instrucao 3:
iexec=1; iprtin=0; irank=0; nsd=2; numnp=$(( (nelDir1+1) * (nelDir2+1) )); ndof=1; nlvect=0; numeg=1; nedge=$(( nelDir1 * (nelDir2+1) + (nelDir1+1) * nelDir2  )); npar=2
printf "%10s %9s %9s %9s %9s %9s %9s %9s %9s %9s \n" ${iexec} ${iprtin} ${irank} ${nsd} ${numnp} ${ndof} ${nlvect} ${numeg} ${nedge} ${npar}

# Instrucao 4 (leitura de coordenadas):
printf "%10s %9s %9s %9s\n" 1 4   ${domI_X} ${domI_Y} 
printf "%10s %9s %9s %9s\n" "" "" ${domF_X} ${domI_Y} 
printf "%10s %9s %9s %9s\n" "" "" ${domF_X} ${domF_Y} 
printf "%10s %9s %9s %9s\n" "" "" ${domI_X} ${domF_Y} 
printf "%10s %9s %9s %9s\n" ${nelDir1} "1" ${nelDir2} $(( nelDir1+1 ))
printf "%10s \n" 0 

# Instrucao 5 (Leitura da localizacao das condicoes de contorno nas arestas):
nleft=$(( nelDir1 * nelDir2 + (nelDir1+1) * (nelDir2-1) + 1 ))
nright=$(( nelDir1 * nelDir2 + (nelDir1+1) * (nelDir2-1) + nelDir1 + 1 ))
incly=$(( 2 * nelDir1 + 1 ))
printf "%10s %9s %9s %9s %9s\n" "1" ${nelDir1} "1" "1" "0"
printf "%10s %9s %9s %9s %9s\n" $((nelDir1 + 1)) ${nleft} ${incly} "1" "0"
printf "%10s %9s %9s %9s %9s\n" $((2 * nelDir1 + 1)) ${nright} ${incly} "1" "0"
printf "%10s %9s %9s %9s %9s\n" $((nedge - nelDir1 + 1)) ${nedge} "1" "1" "0"
#printf "%10s \n" 0
printf "\n"

# Instrucao 6 (Leitura dos dados dos elementos):
ntype=1; numel=$((Elements)); numat=1; nint=16; nen=4; nenp=9; unknown=$((nenp)); npars=2; nints=2;
printf "%10s %9s %9s %9s %9s %9s %9s %9s %9s \n" ${ntype} ${numel} ${numat} ${nint} ${nen} ${unknown} ${nenp} ${npars} ${nints}

# Instrucao 7 (parametros do LDGC):
m=1; delta1=" 2.0d+01"; delta2=" -1.d+00"; gamma=" 0.d+00"; del4="0.00"; theta="-1.0"
printf "%10s %9s %9s %9s %9s %9s \n" ${m} ${delta1} ${delta2} ${gamma} ${del4} ${theta}

# Instrucao 8 (termos gravitacionais):
grav1="1.0"; grav2="0.0"; grav3="0.0"
printf "%10s %9s %9s \n" ${grav1} ${grav2} ${grav3}

# Instrucao 9 (leitura dos valores de contorno prescritos):
#printf " %9s \n" 0

# Instrucao 10 (leitura das conectividades na definicao dos elementos):
n=1; m=1; ng=1
no1=1
no2=$(( no1 + 1 ))
no3=$(( nelDir1+1+1+1    ))
no4=$(( no3-1    ))
printf "%10s %9s %9s %9s %9s %9s %9s\n" ${n} ${m} ${no1} ${no2} ${no3} ${no4} ${ng}
incElemX=1
incNoX=1
incElemY=${nelDir1}
incNoY=$(( nelDir1+1 ))
printf "%10s %9s %9s %9s %9s %9s\n" ${nelDir1} ${incElemX} ${incNoX} ${nelDir2} ${incElemY} ${incNoY}
#printf "%10s \n" 0 
printf "\n" 

# Instrucao 11 (leitura das conectividades para os lados):
n1=1
n2=$((nelDir1+2))
n3=$((2*nelDir1 + 2))
n4=$((nelDir1+1))
printf "%10s %9s %9s %9s %9s %9s\n" ${n} ${n1} ${n2} ${n3} ${n4} ${ng}
printf "%10s %9s %9s %9s %9s %9s\n" ${nelDir1} ${incElemX} ${incNoX} ${nelDir2} ${incElemY} ${incly}

# Instrucao 12 (encerramento):
printf "\n"
printf "%2s %3s \n" 1  0 
printf "*end \n"
