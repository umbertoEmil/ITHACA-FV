from scipy.io import mmread
import matplotlib.pyplot as plt

eigenT = mmread('./ITHACAoutput/POD/Eigenvalues_T')
eigenDeltaT = mmread('./ITHACAoutput/POD/Eigenvalues_deltaT')
eigenLambda = mmread('./ITHACAoutput/POD/Eigenvalues_lambda')

cumEigenT = mmread('./ITHACAoutput/POD/CumEigenvalues_T')
cumEigenDeltaT = mmread('./ITHACAoutput/POD/CumEigenvalues_deltaT')
cumEigenLambda = mmread('./ITHACAoutput/POD/CumEigenvalues_lambda')

plt.figure(1)
plt.xlabel("N of basis")
plt.ylabel("Eigenvalues")
plt.semilogy(eigenT,label='T')
plt.semilogy(eigenDeltaT, label='deltaT')
plt.semilogy(eigenLambda, label='lambda')
plt.legend(loc='upper right')

plt.figure(2)
plt.xlabel("N of basis")
plt.ylabel("Comulative eigenvalues")
plt.plot(cumEigenT,label='T')
plt.plot(cumEigenDeltaT, label='deltaT')
plt.plot(cumEigenLambda, label='lambda')
plt.legend(loc='lower right')

plt.show()
