function a=algebraicConnectivity(adj) 
s=graphSpectrum(adj);
a=s(length(s)-1);
a=a/length(adj);