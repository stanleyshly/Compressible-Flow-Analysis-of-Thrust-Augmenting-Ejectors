'''
purpose : to solve the flow parameters at the end of mixing using the 
 control volume approach in a constant area ejector 

programmer: reimplemented in Python 3 by Stanley Shen

variable key:
	P Static Pressure
	T Temperature
	M Mach Number
	APT Primary to Total Area Ratio
	RHO Density
	PP0 Primary Total Pressure
	PA Secondary Total Temperature 
	TP0 Primary Total Temperature
	TS0 Secondary Stagnation Temperature
	ASP Secondary to Primary Area
	C Speed of Sound
	U Speed of Flow
	MSUB Subsonic Flow Mach Number at the End of Mixing
	MSUP Supersonic Flow Mach Number at the End of Mixing
	M3ID Flow Mach Number for Primary Nozzle
	M3SUB Flow Mach Number at Exit to the Diffuser for Subsonic Solution
	M3SUP Flow Mach Number at Exit to the Diffuser for Supersonic Solution
	TAR Thrust Augmentation Ratio
	subscriptP Primary Flow Conditions at Ejector Inlet
	subscriptS Secondary Flow Conditions at Ejector Inlet 
	subscriptSUB Subsonic Solution to the Mixed Flow
	subscriptSUP Supersonic Solution to the Mixed Flow
	subscript3 Flow conditions at exit to diffuser
'''
import math

K = 1.4
R = 287

#Stagnation Conditions for Primary and Secondary FLows
PP0 = 6# Primary Total Pressure
PA =  1# Secondary Total Pressure
TP0 = 3.35# Primary Total Temperature
TS0 = 1# Secondary Stagnation Temperature
ASP = 10 # Secondary to Primary Area 

# Compute Ideal Flow Mach Number for Primary Nozzle
M3ID=((PP0**((K-1)/K)-1)*2/(K-1))**0.5
T3ID = TP0/(1+((K-1)*M3ID**2)/2)

# Vary Stagnation Inlet pressure from 0.2 atm to 1 atm
PS = 0.2

while (PS < 0.99):
	# Calculate Secondary Flow Mach Number
	MS = (((PA/PS)**((K-1)/K)-1)*2/(K-1))**0.5
	# Calculate Primary Flow Mach Number
	MP = (((PP0/PS)**((K-1)/K)-1)*2/(K-1))**0.5
	
	# Calculate Primary and Secondary Flow Conditions
	TS = TS0/(1+((K-1)/2)*MS**2)
	TP = TP0/(1+((K-1)/2)*MP**2)
	APT = (1/MP**2)*((2/(K+1)) * (1+((K-1)/2)*MP**2))**((K+1)/(K-1))
	
	CS = (K*R*TS)**0.5
	CP = (K*R*TP)**0.5
	US = MS*CS
	UP = MP*CP
	RHOS = (PS/(R*TS))*101300
	RHOP = (PS/(R*TP))*101300
	# Determine Mass Flow Rate Ratio
	MFR = ASP*(MS/MP)*((TP/TS)**0.5)
	# Calculate Pressure and Temperature
	PRP = PS/PP0
	PRS = PS/PA
	TRP = TP/TP0
	TRS = TS/TS0
	N = K/(K-1)
	N1 = (K-1)/K
	N2 = (K-1)/2

	J = ((TRP **0.5)*(((ASP+1) / (K*MP) ) +MP) + MFR*MS*((TS/TP0)**0.5) \
		/((1+MFR)*(TS0/TP0)) * (1+MFR))**0.5
	A = 1-(J**2)*(K-1)/2
	B = 2-K*(J**2)
	DET = (B**2)-(4*A)
	if (DET < 0):
		print('imaginary solution')
	else:
		#compute solutions to the flow at the end of mixing
		MSUB = ((-B-(DET**0.5))/(2*K*A))**0.5
		MSUP = ((-B+(DET**0.5))/(2*K*A))**0.5
		#calculate flow conditions at the end of mixing on the subsonic branch
		TSUB0 = (TP0 + MFR*TS0)/(1+MFR)
		XSUB = (1+(N2*(MP**2)))/(1+(N2*(MSUB**2)))
		ZSUB = (1+MFR)*MP*((TSUB0*XSUB/TP0)**0.5)/((ASP+1)*MSUB)
		TSUB = TS*XSUB*TSUB0/TS0
		PSUB = PS*ZSUB
		PSUB0 = ((1+N2*(MSUB*2))**N)*PSUB
		YSUB = TSUB/TP
		WSUB = TSUB/TS0
		Y1 = YSUB
		W1 = WSUB
		Z1 = ZSUB
		DSSUB = N * math.log(Y1) + N*MFR*math.log(W1)-(1+MFR)*math.log(Z1) 
		# Calculate Diffuser Exit Flow Conditions for Subsonic Mixed Flow
		N3 = PSUB0**N1
		N4 = 2*((PSUB0**N1)-1)
		M3SUB = (2*((PSUB0**N1)-1)/(K-1))**0.5
		T3SUB = TSUB0/(1+((K-1)*M3SUB**2)/2)
		# Evaluate Thrust Augmentation Ratio for Subsonic Mixed Flow
		TARSUB = (1+MFR)*M3SUB*(T3SUB**0.5)/(M3ID*(T3ID**0.5)) # possible error here
		# Calculate Flow Conditions for Supersonic Mixed Flow at End of Mixing
		XSUP = (1+N2*(MP**2))/(1+N2*(MSUP**2))
		ZSUP = (1+MFR) * MP * ((TSUB*XSUP/TP0)**0.5)/((ASP+1)*MSUP)
		PSUP = PS*ZSUP
		PSUP0 = ((1+((K-1)*MSUP**2)/2)**N)*PSUP
		TSUP = TS*XSUP*(TSUB0/TS0)
		YSUP = TSUP/TP
		WSUP = TSUP/TS
		#print(YSUP, WSUP)
		try:
			# Calculate Flow Conditions for Supersonic Mixed Flow At Exit to the Diffuser

			DSSUP = N* math.log(YSUP)+N*MFR*math.log(WSUP)-(1+MFR)*math.log(ZSUP) # possible error here
			M3SUP = (2*((PSUP0**(1/N))-1)/(K-1))**0.5
			T3SUP = TSUB0/(1+((K-1)*M3SUP**2)/2)
			# Evaluate the Thrust Augmentation Ratio For the Supersonic Solution
			TARSUP = (1+MFR) *M3SUP*(T3SUB**0.5)/(M3ID*(T3ID**0.5))
			print(TARSUB, TARSUP)
		except: 
			#print('imaginary solution')
			print(TARSUB)


	PS = PS + 0.02