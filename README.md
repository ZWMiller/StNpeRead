# StNpeRead
StNpeRead Online Code
Z. Miller - 8/3/15

This takes in the StDmesonTreeMaker output root files (Mustafa Mustafa) as input. The StDmestonTreeMaker selects, with loose cuts: electrons, kaons, pions, electron pairs, etc. These are sorted into their own branches. The important branch for this analysis is the "All tracks" branch which allows all identified electron and hadron tracks into the branch. (NPE = non-photonic electron (non-pair))

Using this input, this code sorts through all tracks looking for electrons, then sorts through all tracks in the same event as an electron looking for hadron type events and calculates delta_Phi between the two. It also does this for all electron pairs (which come from photon conversion, eta-Dalitz and pi-Dalitz), using the primary electron as the "electron phi." The pairs can be separated into like sign pairs and unlike sign pairs.

Final result is: 
delPhi_NPE = delPhi_semiInclusive*purity_sample - 1/epsilon(delPhi_US - delPhi_LS)
epsilon = photonic reconstruction efficiency. US-LS = photonic electrons. Semi-Inclusive is all electrons except reconstructed photonics.
