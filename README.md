# disordered_protein_polymer_chain_theory
in this project, i am attempting to use polymer physics quantities to describe and generate a map of the conformational landscapes of intrinsically disordered proteins. This map shows what the range of conformations accessible to a protein is in terms of the polymer physics quantities of shape and size. 
We calculate a quantity we call the instantaneous shape ratio - this quantity is the square of end-to-end-distance distance divided by the square of radius of gyration. This is our measure of protein shape. 
We use the protein's radius of gyration as a measure of the protein's size.

i have Monte Carlo codes for the Random Walk (RW), Self-Avoiding Walk (SAW) and Gaussian Walk (GW) polymer chain models. i compute shape ratio values by taking square of end-to-end distance and dividing by square of radius of gyration. this shape ratio value is compbined with radius of gyration to generate a map of conformational landscapes. Maps of select proteins are plotted against GW polymer chain model. 
