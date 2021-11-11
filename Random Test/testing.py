from Comparing import compare
mxpera3=-1
mxpera=-1
mnpera3=100000
mnpera=100000
for i in range(100):
    mxper3,mxper,mnper3,mnper=compare()
    mxpera3=max(mxpera3,mxper3)
    mxpera=max(mxpera,mxper)
    mnpera=min(mnper,mnpera)
    mnpera3=min(mnper3,mnpera3)
    print(i)

print("Comparing with PageRank Algorithm")
print("max percentage error", mxpera3," min percentage error ", mnpera3)
print("Comparing with Power Iteration Algorithm")
print("max percentage error ", mxpera," min percentage error ", mnpera)