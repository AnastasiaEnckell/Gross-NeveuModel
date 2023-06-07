# This program simulate a real-time dynimics of the massive Gross-Neveu model with the Bravyi-Kitaev transformation and a lexicographic ordering of terms. 

using CSV,  DataFrames
include("QSim.jl")


N = 16+1 # total number of qubits
maxdim = 5 # maximum allowed bond dimension for MPS initiallly 20
dt = 0.1 #timestep 
lastqubit = N

#some constants
r = 1
k = 1
m = 1
g = 1

backend_MPS = "MPS_ITensor" # calculate with MPS (any higher number of qubits)


qc_MPS = initialise_qcircuit(N, backend_MPS, maxdim=maxdim)


df = CSV.read("xyzdata.csv", DataFrame; header=false) #open the file with the data 
vectorOfVectors = [collect(row) for row in eachrow(df)]
new_length = length(vectorOfVectors) ÷ 3

#Reshape the data into the suitable form 
u = [vectorOfVectors[(i-1)*3+1:i*3] for i in 1:new_length] 
gates = [[[x for x in inner_vec if x != 0] for inner_vec in mid_vec] for mid_vec in u]


dg = CSV.read("zdata.csv", DataFrame; header=false)
vectorOfVectors = [collect(row) for row in eachrow(dg)]
self = [[x for x in inner_vec if x != 0] for inner_vec in vectorOfVectors]


#Example of an initial particle on a site 8, calculated with the help on the Bravyi-Kitaev matrix. 
PauliX!(qc_MPS, [8])
PauliX!(qc_MPS, [16])

for t = 0:3 #timesteps
  for i=1:length(gates)
        xyz = gates[i]

        if i > 1
          xyzbefore = gates[i-1]
        else
          xyzbefore = [0,0,0]
        end 

        if i < length(gates)
          xyzafter = gates[i+1]
        else 
          xyzafter = [0,0,0]
        end


        for element in sort(xyz[1])
          if element ∉ xyzbefore[1]
              Hadamard!(qc_MPS, [element])
           end
        end

        for element in sort(xyz[2])
      α = Complex(1)
      β = 0
      γ = 0
      δ = -im
            if element ∉ xyzbefore[2]
              apply_single_site_gates!(qc_MPS, [element], "U", α, β, γ, δ)
              Hadamard!(qc_MPS, [element])
           end
        end
        
        for element in sort(xyz[1])
            if element ∉ xyzbefore[1]
                Cnot!(qc_MPS, [element, lastqubit])
            end
        end

        for element in sort(xyz[2])
            if element ∉  xyzbefore[2]
                Cnot!(qc_MPS, [element, lastqubit])
            end
        end

        for element in sort(xyz[3])
           if element > 0 && element ∉  xyzbefore[3]
                Cnot!(qc_MPS, [element, lastqubit])
           end
        end


      if -99 in xyz[3]
        if -1 in xyz[3]
            θ = -1/2*dt
            RzGate!(qc_MPS, [lastqubit], θ)
          end
          if -2 in xyz[3]
            θ = 1/2*dt
            RzGate!(qc_MPS, [lastqubit], θ)
          end
          if -3 in xyz[3]
            θ = -r/2*dt
            RzGate!(qc_MPS, [lastqubit], θ)
          end
          if -4 in xyz[3]
            θ = r/2*dt
            RzGate!(qc_MPS, [lastqubit], θ)
          end
        else
            if -1 in xyz[3]
              θ = 1/2*dt
              RzGate!(qc_MPS, [lastqubit], θ)
            end
            if -2 in xyz[3]
              θ = -1/2*dt
              RzGate!(qc_MPS, [lastqubit], θ)
            end
            if -3 in xyz[3]
              θ = r/2*dt
              RzGate!(qc_MPS, [lastqubit], θ)
            end
            if -4 in xyz[3]
              θ = -r/2*dt
              RzGate!(qc_MPS, [lastqubit], θ)
            end
        end
        
        for element in reverse(sort(xyz[3]))
            if element > 0 && element ∉  xyzafter[3]
                Cnot!(qc_MPS, [element, lastqubit])
            end
        end

        for element in reverse(sort(xyz[2]))
            if element ∉  xyzafter[2]
                Cnot!(qc_MPS, [element, lastqubit])
            end
        end

        for element in reverse(sort(xyz[1]))
            if element ∉  xyzafter[1]
                Cnot!(qc_MPS, [element, lastqubit])
            end
        end

        for element in reverse(sort(xyz[2]))
      α = Complex(1)
      β = 0
      γ = 0
      δ = im
          if element ∉  xyzafter[2]
              Hadamard!(qc_MPS, [element])
              apply_single_site_gates!(qc_MPS, [element], "U", α, β, γ, δ)
          end
        end

        for element in reverse(sort(xyz[1]))
          if element ∉  xyzafter[1]
              Hadamard!(qc_MPS, [element])
            end
        end

    end


  for j = 1:length(self)
      
      if j > 1
          xyzbefore = self[j-1]
        else
          xyzbefore = [0,0,0]
        end 

        if j < length(self)
          xyzafter = self[j+1]
        else 
          xyzafter = [0,0,0]
        end


    for index in self[j]
      if index > 0
        if index ∉ xyzbefore
            Cnot!(qc_MPS, [index, lastqubit])
              end
            end
      end

      if -1 in self[j]
        θ = g*k*dt/2
        RzGate!(qc_MPS, [lastqubit], θ)
      end


      if -2 in self[j]
        θ = -g*k*dt/2
        RzGate!(qc_MPS, [lastqubit], θ)
      end
      if -3 in self[j]
          θ = (m*k+r)*dt + (g)*k*dt/2
        RzGate!(qc_MPS, [lastqubit], θ)
      end
      if -4 in self[j]
          θ = -(m*k+r)*dt + (g)*k*dt/2
        RzGate!(qc_MPS, [lastqubit], θ) 
      end

    for index in reverse(self[j])
      if index > 0
        if index ∉ xyzafter
            Cnot!(qc_MPS, [index, lastqubit])
              end
            end
      end

  end


  N_meas = 100  
  register1 = [i for i in 1:N-1]
  meas1 = sample_measurement(qc_MPS, register1, N_meas, save_measurement=false)


  open("bravyi-kitaev.txt", "a") do io
      writedlm(io, meas1)
      writedlm(io, " ")
  end

end

draw(qc_MPS)

