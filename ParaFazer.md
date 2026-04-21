# Lista de funcionalidades para implementação

## Sistemas Exemplo
- [ ] Implementação do sistema massa-mola.
  - [ ] Adição da possibilidade de animação.
  - [ ] Adição dos termos não lineares.
  - [ ] Adição da possibilidade de controle.
  
## FrameWork
- [ ] Melhorar a documentação do framework
- [ ] Melhor separação das funções de controle das plantas.
- [ ] Adição de funções para aplicação de controladores em série com a planta.
- [ ] Adição de funções para aplicação de entrada/disturbio no sistema
  
### Ferramentas
- [ ] Adição do de processamento paralelo na CPU e GPU.
  - [x] Criação de templates para uso da GPU.
  - [x] Criação de templates para uso dos núcleos da CPU.
  - [ ] Implementação de funções que usam esses recursos.
- [ ] Adição de visualização do retrato de fase de um sistema.
  - [x] Possibilidade de visualização de uma solução única.
  - [ ] Possibilidade de visualização de todo o retrado para um dado sistema.
  
### Controladores
- [ ] Adição do controlador PID.
  - [ ] Deixar o PID de uma forma mais simples de ser usada.
- [ ] Adição do controlador LQR. 
- [ ] Adição do controlador MPC.
- [ ] Adição do controlador Fuzzy. 

## Drone - TCC
- [ ] Adicionar um editor latex 
- [ ] Pensar em uma estrutura padrão para aplicação da física e dos controladores diferentes.
- [ ] Modificar a Iteração DK para conseguir implementar o ponto abaixo.
- [ ] Calcular a função de controle de acordo com o ponto de operação do drone, colocar iteracaoDK dentro da função de controle.
- [ ] Conseguir provar que o controlador funciona para uma região de operação.
- [ ] Adição de um algoritmo de otimização para os controladores achados.
- [ ] Implementar seguimento de trajetória e aplicar o cálculo do novo controlador de acordo com os novos pontos de operação.

