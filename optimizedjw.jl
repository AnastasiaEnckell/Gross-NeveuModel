# This program simulate a real-time dynimics of the massive Gross-Neveu model with the Jordan-Wigner transformation and a lexicographic ordering of terms. 

include("QSim.jl")


N =16+1 # total number of qubits
leng = round(Int, (N-1)/4 - 1)
maxdim = 3 # maximum allowed bond dimension for MPS initiallly 20
dt = 0.1 #timestep 
last_qubit = N

#some constants
r = 1
a = 1
m = 1
g = 1

backend_MPS = "MPS_ITensor" # calculate with MPS (any higher number of qubits)

qc_MPS = initialise_qcircuit(N, backend_MPS, maxdim=maxdim)


PauliX!(qc_MPS, [8]) #an initial particle at the on the position 8


#quantum circuit 
for t = 0:10 #timesteps

	for j in 0:leng

		α = Complex(1)
		β = 0
		γ = 0
		δ = -im


		α1 = Complex(1)
		β1 = 0
		γ1 = 0
		δ1 = im



		######1st term

		Hadamard!(qc_MPS, [1 + 4j, 1 + mod(4+4j, N-1)])


		if j < 1
			for i in 0:4
				Cnot!(qc_MPS, [1 + mod(i + 4j, N-1),last_qubit])
			end
		else 
			Cnot!(qc_MPS, [1 + mod(0 + 4j, N-1),last_qubit])
			Cnot!(qc_MPS, [1 + mod(3 + 4j, N-1),last_qubit])
			Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1),last_qubit])
		end 

		θ = r/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)



		Cnot!(qc_MPS, [ 1 + mod(4 + 4j, N-1), last_qubit])


		Hadamard!(qc_MPS, [1 + 4j, 1 + mod(4+4j, N-1)])





		######2nd term

		Hadamard!(qc_MPS, [1+4j, 3 + mod(4+4j, N-1)])



		Cnot!(qc_MPS, [1+ mod(4 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1+ mod(5 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1+ mod(6 + 4j, N-1), last_qubit])


		θ = -1/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)


		Cnot!(qc_MPS, [1+ mod(6 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1+ mod(5 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1+ mod(4 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1+ mod(0 + 4j, N-1),last_qubit])

		
		Hadamard!(qc_MPS, [1+4j, 3 + mod(4+4j, N-1)])






		########3rd term

		apply_single_site_gates!(qc_MPS, [1 + 4j, 1 + mod(4+4j, N-1)], "U", α, β, γ, δ)
		Hadamard!(qc_MPS, [1 + 4j, 1 + mod(4+4j, N-1)])


		Cnot!(qc_MPS, [1 + mod(0 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1),last_qubit])


		θ = r/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)


		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1),last_qubit])

		Hadamard!(qc_MPS, [1 + 4j, 1 + mod(4+4j, N-1)])
		apply_single_site_gates!(qc_MPS, [1 + 4j, 1 + mod(4+4j, N-1)], "U", α1, β1, γ1, δ1)



		#########4th term


		apply_single_site_gates!(qc_MPS, [1+4j, 3 + mod(4+4j, N-1)], "U", α, β, γ, δ)
		Hadamard!(qc_MPS, [1+4j, 3 + mod(4+4j, N-1)])


		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1),last_qubit ])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit ])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1),last_qubit ])
	

		θ = -1/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)


		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(1 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(0 + 4j, N-1), last_qubit])


		Hadamard!(qc_MPS, [1+4j, 3 + mod(4+4j, N-1)])
		apply_single_site_gates!(qc_MPS, [1+4j, 3 + mod(4+4j, N-1)], "U", α1, β1, γ1, δ1)








		######5th term

		Hadamard!(qc_MPS, [2 + 4j, 2 + mod(4+4j, N-1)])


		Cnot!(qc_MPS, [1 + mod(1 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1+ mod(5 + 4j, N-1),last_qubit])

		θ = r/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit])


		Hadamard!(qc_MPS, [2 + 4j, 2 + mod(4+4j, N-1)])




		#########6th term

		
		Hadamard!(qc_MPS, [2+4j, 4 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(7 + 4j, N-1),last_qubit])

		θ = -1/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)


		Cnot!(qc_MPS, [1 + mod(7 + 4j, N-1), last_qubit ])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1), last_qubit ])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit ])
		Cnot!(qc_MPS, [1 + mod(1 + 4j, N-1), last_qubit ])

		Hadamard!(qc_MPS, [2+4j, 4 + mod(4+4j, N-1)])





		########7th term

		apply_single_site_gates!(qc_MPS, [2 + 4j, 2 + mod(4+4j, N-1)], "U", α, β, γ, δ)
		Hadamard!(qc_MPS, [2 + 4j, 2 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(1 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit])

		θ = r/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])

		Hadamard!(qc_MPS, [2 + 4j, 2 + mod(4+4j, N-1)])
		apply_single_site_gates!(qc_MPS, [2 + 4j, 2 + mod(4+4j, N-1)], "U", α1, β1, γ1, δ1)




		#########8th term

		apply_single_site_gates!(qc_MPS, [2+4j, 4 + mod(4+4j, N-1)], "U", α, β, γ, δ)
		Hadamard!(qc_MPS, [2+4j, 4 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(7 + 4j, N-1),last_qubit])

		θ = -1/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(1 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(2 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(7 + 4j, N-1),last_qubit])

		Hadamard!(qc_MPS, [2+4j, 4 + mod(4+4j, N-1)])
		apply_single_site_gates!(qc_MPS, [2+4j, 4 + mod(4+4j, N-1)], "U", α1, β1, γ1, δ1)





		#######9th term 

		Hadamard!(qc_MPS, [3+4j, 1 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(2 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1),last_qubit])

		θ = 1/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1),last_qubit])
		
		Hadamard!(qc_MPS, [3+4j, 1 + mod(4+4j, N-1)])






		#######10th term 

		Hadamard!(qc_MPS, [3+4j, 3 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1), last_qubit])

		θ = -r/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(2 + 4j, N-1), last_qubit])

		Hadamard!(qc_MPS, [3+4j, 3 + mod(4+4j, N-1)])



		########11th term 

		apply_single_site_gates!(qc_MPS, [3+4j, 1 + mod(4+4j, N-1)], "U", α, β, γ, δ)
		Hadamard!(qc_MPS, [3+4j, 1 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(2 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1), last_qubit])

		θ = 1/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [ 1 + mod(4 + 4j, N-1),last_qubit])

		Hadamard!(qc_MPS, [3+4j, 1 + mod(4+4j, N-1)])
		apply_single_site_gates!(qc_MPS, [3+4j, 1 + mod(4+4j, N-1)], "U", α1, β1, γ1, δ1)




		########12th term

		apply_single_site_gates!(qc_MPS, [3+4j, 3 + mod(4+4j, N-1)], "U", α, β, γ, δ)
		Hadamard!(qc_MPS, [3+4j, 3 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1), last_qubit])

		θ = -r/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(3 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(2 + 4j, N-1), last_qubit])

		Hadamard!(qc_MPS, [3+4j, 3 + mod(4+4j, N-1)])
		apply_single_site_gates!(qc_MPS, [3+4j, 3 + mod(4+4j, N-1)], "U", α1, β1, γ1, δ1)




		######13th term

		Hadamard!(qc_MPS, [4+4j, 2 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(3 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])

		θ = 1/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit])

		Hadamard!(qc_MPS, [4+4j, 2 + mod(4+4j, N-1)])
		#count = count + 2






		#######14th term

		Hadamard!(qc_MPS, [4+4j, 4 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(7 + 4j, N-1), last_qubit])

		θ = -r/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(7 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(3 + 4j, N-1), last_qubit])

		Hadamard!(qc_MPS, [4+4j, 4 + mod(4+4j, N-1)])





		#######15th term

		apply_single_site_gates!(qc_MPS, [4+4j, 2 + mod(4+4j, N-1)], "U", α, β, γ, δ)
		Hadamard!(qc_MPS, [4+4j, 2 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(3 + 4j, N-1),last_qubit])
		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])

		θ = 1/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1),last_qubit])

		Hadamard!(qc_MPS, [4+4j, 2 + mod(4+4j, N-1)])
		apply_single_site_gates!(qc_MPS, [4+4j, 2 + mod(4+4j, N-1)], "U", α1, β1, γ1, δ1)

	


		
		########16th term

		apply_single_site_gates!(qc_MPS, [4+4j, 4 + mod(4+4j, N-1)], "U", α, β, γ, δ)
		Hadamard!(qc_MPS, [4+4j, 4 + mod(4+4j, N-1)])

		Cnot!(qc_MPS, [1 + mod(5 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(6 + 4j, N-1), last_qubit])
		Cnot!(qc_MPS, [1 + mod(7 + 4j, N-1), last_qubit])

		θ = -r/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)

		if j == leng
		for i in reverse(3:7)
			Cnot!(qc_MPS, [1 + mod(i + 4j, N-1), last_qubit])
			count = count +1
		end
		else 
			Cnot!(qc_MPS, [1 + mod(3 + 4j, N-1), last_qubit])
			Cnot!(qc_MPS, [1 + mod(4 + 4j, N-1), last_qubit])
			Cnot!(qc_MPS, [1 + mod(7 + 4j, N-1), last_qubit])
		end


		Hadamard!(qc_MPS, [4+4j, 4 + mod(4+4j, N-1)])
		apply_single_site_gates!(qc_MPS, [4+4j, 4 + mod(4+4j, N-1)], "U", α1, β1, γ1, δ1)


end



for j in 0:leng
		#######17

		Cnot!(qc_MPS, [1+4j,last_qubit])
		Cnot!(qc_MPS, [2+4j,last_qubit])
		θ = -a*g/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [2+4j,last_qubit])

		Cnot!(qc_MPS, [ 3+4j,last_qubit])
		θ = a*g/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [ 3+4j,last_qubit])

		Cnot!(qc_MPS, [ 4+4j,last_qubit])
		θ = a*g/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [ 4+4j,last_qubit])
		Cnot!(qc_MPS, [ 1+4j,last_qubit])

		Cnot!(qc_MPS, [ 2+4j, last_qubit])
		Cnot!(qc_MPS, [3+4j, last_qubit])
		θ = a*g/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [3+4j, last_qubit])		

		Cnot!(qc_MPS, [4+4j, last_qubit])
		θ = a*g/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [ 2+4j, last_qubit])

		Cnot!(qc_MPS, [ 3+4j, last_qubit])
		θ = -a*g/2*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [ 4+4j, last_qubit])
		Cnot!(qc_MPS, [ 3+4j, last_qubit])

		Cnot!(qc_MPS, [ 1+4j, last_qubit])
		θ = (m*a+r)*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [1+4j, last_qubit])

		Cnot!(qc_MPS, [2+4j, last_qubit])
		θ = (m*a+r)*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [2+4j, last_qubit])

		Cnot!(qc_MPS, [3+4j, last_qubit])
		θ = -(m*a+r)*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [3+4j, last_qubit])

		Cnot!(qc_MPS, [4+4j, last_qubit])
		θ = -(m*a+r)*dt
		RzGate!(qc_MPS, [last_qubit], θ)
		Cnot!(qc_MPS, [4+4j, last_qubit])

	end

	
	N_meas = 100

	register1 = [i for i in 1:N-1]
	meas1 = sample_measurement(qc_MPS, register1, N_meas, save_measurement=false)

	open("32bd15kr0.001p1418.txt", "a") do io
           writedlm(io, meas1)
           writedlm(io, " ")
     end
	
end


draw(qc_MPS)













