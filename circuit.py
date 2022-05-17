from PySpice.Spice.NgSpice import Shared
import PySpice.Unit as PU
from PySpice.Spice.Netlist import Circuit
import numpy as np
import random
import copy


elements = {  "R1": (1, 500@PU.u_Ω),
              "R2": (2, 10@PU.u_kΩ),
              "R3": (3, 1@PU.u_kΩ),
              "R4": (4, 20@PU.u_kΩ),
              "C1": (1, 1@PU.u_uF),
              "C2": (2, 5@PU.u_uF),
              "C3": (3, 7@PU.u_uF),
              "L1": (1, 1@PU.u_uH),
              "L2": (2, 10@PU.u_uH),
              "L3": (3, 15@PU.u_uH)
              }

class Circus:

    def __init__(self, elements_list, N_orig_comp,gene=None):
        self.elements_list = elements_list
        self.elements_names = list(elements_list.keys())
        self.N_orig_comp = N_orig_comp
        self.create_node_list()
        if gene==None:
            self.generate_gene()
        else:
            self.gene=gene
            self.list_nodes_comps()
        # self.add_jumpers()
        # self.draw_circuit()

    def create_node_list(self):
        
        self.node_list = ["Vin"]
        
        for i in range(1, self.N_orig_comp+1):
            n_node_list = f"n{i}"
            self.node_list.append(n_node_list)
            
        self.node_list.append("GND")
        

    def generate_gene(self):
        # without Vin and GND:
        self.out = random.choice(self.node_list[1:-1])

        # with Vin and GND:
        # self.out = random.choice(self.node_list)

        self.components = random.choices(self.elements_names, k=self.N_orig_comp)
        self.nodes = []
        for comp in self.components:
            start = np.random.randint(0, self.N_orig_comp)
            end = np.random.randint(start + 1, self.N_orig_comp + 2)
            self.nodes.append(self.node_list[start])
            self.nodes.append(self.node_list[end])
        self.assemble_gene()

    def assemble_gene(self):
        self.gene = [self.out]
        for c, comp in enumerate(self.components):
            self.gene.append(self.nodes[2 * c])
            self.gene.append(comp)
            self.gene.append(self.nodes[2 * c+1])

    def list_nodes_comps(self):
        self.nodes=[]
        self.components=[]
        for i in range(self.N_orig_comp):
            self.nodes.append(self.gene[i*3+1])
            self.components.append(self.gene[i*3+2])
            self.nodes.append(self.gene[i*3+3])

    def add_jumpers(self):
        #IMPORTANT COMMENT: this function add jumpers
        # check which node is first and last, and change Vin and GND accordingly
        # exclude out position from the search
        self.node_list_new=copy.copy(self.node_list)
        primeiro="Vin"
        for node in self.node_list:
            if node in self.gene[1:]:
                primeiro=node
                break
        ultimo="GND"
        for node in reversed(self.node_list):
            if node in self.gene[1:]:
                ultimo=node
                break
        print(primeiro,ultimo)
        self.gene_jumper=["Vin" if node==primeiro else node for node in self.gene]
        self.gene_jumper=["GND" if node==ultimo else node for node in self.gene_jumper]
        
        #check for empty nodes, and delete them:
        for node in self.node_list:
            if node not in self.gene_jumper:
                self.node_list_new.remove(node)

        # OUT position. If at Vin or GND, it changes places:
        troca=self.gene_jumper[0]
        if troca=="Vin":
            self.gene_jumper[0]=self.node_list_new[1]
        elif troca=="GND":
            self.gene_jumper[0]=self.node_list_new[-2]

        #check paths:
        # gene=self.gene_jumper[1:self.N_orig_comp*3+1]
        # indices = [i for i, x in enumerate(gene) if x == "Vin"]
        # for ind in indices:
        #     print(gene[ind+2])
        #     GND ENCERRA O WHILE. FAZER WHILE ATÉ TERMINAR O CAMINHO
        #     if gene[ind+2]=="GND":
        #         ignore=True
        #     else:
        #         indices2 = [i for i, x in enumerate(gene) if x == gene[ind+2]]
        #         for ind2 in indices2:
        #             print(ind2)

        # add jumpers in a mysterious way
        conns_circuit={}
        for i in range(self.N_orig_comp):
            node1=self.gene_jumper[i*3+1]
            component=self.gene_jumper[i*3+2]
            node2=self.gene_jumper[i*3+3]
            node=[node1,node2]
            node.sort()
            if (node1+node2 in conns_circuit):
                conns_circuit[node1+node2].append(component)
            else:
                conns_circuit[node1+node2]=[component]
        
        for n in range(len(self.node_list_new)-1):
            node1 = self.node_list_new[n]
            node2 = self.node_list_new[n+1]
            if node1+node2 not in conns_circuit:
                self.gene_jumper.append(node1)
                self.gene_jumper.append("jumper")
                self.gene_jumper.append(node2)


    def assemble_circuit(self,verbose=False):
        self.circuit = Circuit("pindamonhangaba")
        self.circuit.SinusoidalVoltageSource('input', 'Vin', self.circuit.gnd, amplitude=1@PU.u_kV)
        K=len(self.gene_jumper)
        # OUT position:
        troca=self.gene_jumper[0]

        new_gene=["out" if comp==troca else comp for comp in self.gene_jumper]
        for pos in range(1,K,3):
            # print(pos,self.gene[pos:pos+3])
            pos1=new_gene[pos]
            pos2=new_gene[pos+2]
            if pos1!=pos2:
                component=new_gene[pos+1]
                if component=="jumper":
                    func_add=self.circuit.R
                    valor=1@PU.u_mΩ
                else:
                    if (component[0] == "R"):
                        func_add=self.circuit.R
                    elif (component[0] == "C"):
                        func_add=self.circuit.C
                    elif (component[0] == "L"):
                        func_add=self.circuit.L
                    valor=self.elements_list[component][1]
                func_add(pos,pos1,pos2,valor)

    def draw_circuit(self,draw_jumpers=True):
        self.node_list_dict={}
        for n,node in enumerate(self.node_list_new):
            self.node_list_dict[node]=n

        conns={}
        print("    ".join(self.node_list_new))
        pos=self.node_list_dict[self.gene_jumper[0]]
        print(" "+" "*pos*6+"*")
        circuit=self.gene_jumper[1:]
        if draw_jumpers:
            N=int((len(self.gene_jumper)-1)/3)
        else:
            N=self.N_orig_comp
        refugo=[]
        for i in range(N):
            node1=circuit[i*3]
            component=circuit[i*3+1].replace("jumper","jp")
            node2=circuit[i*3+2]
            pos1=self.node_list_dict[node1]
            pos2=self.node_list_dict[node2]
            #check if pos1 and pos2 are inverted:
            if pos1>pos2:
                pos_inter=pos1
                pos1=pos2
                pos2=pos_inter
            dist=((pos2-pos1)*6)-2
            dist2=int(dist/2)
            if pos1!=pos2:
                print("  "+pos1*6*" "+"|"+(dist2-1)*"-"+component+(dist2-1)*"-"+"|")
            else:
                # print(pos1,pos2,component,node1,node2)
                refugo.append(component)
        if len(refugo)!=0:
            print("ligados em si mesmo:", refugo)


    def mutation(self):
        porcentagem=5
        for c in range(self.N_orig_comp):
            aux = random.randint(0, 100)
            if aux <= porcentagem:
                self.components[c] = random.choice(self.elements_names)
        for n in range(self.N_orig_comp * 2):
            aux = random.randint(0, 100)
            if aux <= porcentagem:
                self.nodes[n] = random.choice(self.node_list)
        self.assemble_gene()


def single_cross(circ1, circ2, verbose=False):
    N = circ1.N_orig_comp
    cut = random.randint(0, len(circ1.gene))
    cross1 = circ1.gene[:cut] + circ2.gene[cut:]
    cross2 = circ2.gene[:cut] + circ1.gene[cut:]

    if verbose:
        print(circ1.gene, circ2.gene)
        print(cut)
        print(cross1, cross2)

    return cross1, cross2


def double_cross(circ1, circ2, verbose=False):
    N = circ1.N_orig_comp
    cut1 = random.randint(0, len(circ1.gene) - 2)
    cut2 = random.randint(cut1 + 1, len(circ1.gene))

    cross1 = circ1.gene[:cut1] + circ2.gene[cut1:cut2] + circ1.gene[cut2:]
    cross2 = circ2.gene[:cut1] + circ1.gene[cut1:cut2] + circ2.gene[cut2:]

    if verbose:
        print(circ1.gene, circ2.gene)
        print(cut1, cut2)
        print(cross1, cross2)

    return cross1, cross2


def crossover(circ1, circ2, verbose):
    chance = np.random.randint(0, 100)
    if chance < 30:
        cross1, cross2 = single_cross(circ1, circ2, verbose)
    else:
        cross1, cross2 = double_cross(circ1, circ2, verbose)
    new_circ1=Circus(circ1.elements_list,circ1.N_orig_comp,cross1)
    new_circ2=Circus(circ2.elements_list,circ2.N_orig_comp,cross2)
    new_circ1.out=circ1.out
    new_circ2.out=circ2.out
    return new_circ1,new_circ2


def eval_fitness(circuit):
    try:
        simulator = circuit.simulator(temperature=25, nominnal_temperature=25)
        analysis = simulator.ac(start_frequency=150@PU.u_Hz, stop_frequency=3@PU.u_kHz, number_of_points=1, variation='dec')
        gain = 20 * np.log10(np.absolute(analysis.out)/(1@PU.u_V))
        delta_G = ((gain[1]) - (gain[0]))
        delta_F = (np.log10(3000) - np.log10(150))
        inclination = delta_G / delta_F
        return inclination
    except NameError:
        return 1e9


if __name__ == "__main__":
    verb = True
    N = 5
    circ1 = Circus(elements, N)
    # circ2=Circus(elements,N)
    # cross1,cross2=crossover(circ1,circ2,verb)
    # print("original gene")
    # print(circ1.gene)
    # print(circ1.nodes)
    # circ1.mutation()
    print(circ1.gene)
    circ1.add_jumpers()
    # circ1.assemble_circuit(verb)
    # fitness=eval_fitness(circ1.circuit)
    # print(fitness)
    print(circ1.gene_jumper)
    circ1.draw_circuit()
