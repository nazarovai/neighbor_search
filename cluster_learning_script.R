# Neighbor search algorithm

#PRELIMINARIES

library(combinat) #for regular expressions

library(mclust) #for mixture of Gaussians cluster analysis

options(stringsAsFactors = F)

#USER-DEFINED VARS

cycles <- 15 #the number of cycles the simulation will go through

phon_dist_threshold <- 1 #the number of traditional features in which a "neighboring" segment may differ from a prototype segment in the neighborhood search algorithm

accuracy_threshold <- 0.01 #this is the max. gain value a constraint must minimally have to be considered "relevant"


#REPRESENTATIONAL SPACE

#define sigma

sigma = c("a","i","u","p","t","k","b","d","g","m","n","q") # q = engma  #The segments allowed in the language.

#Define the space of representations (omega) that the grammar will consider
#This will be all possible CVCVC strings with segments drawn from sigma

all.repr <- expand.grid(							#a data frame of all possible CVCVCs
				c("p","t","k","b","d","g","m","n","q"),
				c("a","i","u"),
				c("p","t","k","b","d","g","m","n","q"),
				c("a","i","u"),
				c("p","t","k","b","d","g","m","n","q")
				)


#Now make the rows of this data frame into strings, and put them all into a vector

omega <- unlist( apply( all.repr, 1, paste, sep="", collapse="" ) ) #contains 248832 representations total


#Now compute the set of all observed forms given the setup of the toy language

all.obs <- expand.grid(
				c("p","t","k","b","d","g"), #first position contains non-nasal consonants
				c("a","i","u"), #second position contains only vowels
				c("p","t","k","b","d","g","m","n","q"), #third position may contain any consonant
				c("a","i","u"), #fourth position contains only vowels
				c("p","t","k","b","d","g","n","q") #fifth position contains all consonants except "m"
					)

all.obs <- subset(all.obs, !(all.obs[,2]=="u" & all.obs[,3] %in% c("p","b","m") & all.obs[,4]=="u") ) #Remove all "u"-labial-"u" sequences from this frame

#Now, put everything in a matrix of strings

observed <- unlist( apply( all.obs, 1, paste, sep="", collapse="" ) )

#PHONETIC DISTANCE

# define phonetic distance between members of sigma (in the form of a data frame, with both row and column names equal to the set of segments in sigma)

phonetic_distance = 	data.frame(
				"a"=c(0,3,3,rep(5,9)),
				"i"=c(3,0,3,rep(5,9)),
				"u"=c(3,3,0,rep(5,9)),

				"p"=c(rep(5,3), 0,1,2, 1,2,3, 2,3,4),
				"t"=c(rep(5,3), 1,0,1, 2,1,2, 3,2,3),
				"k"=c(rep(5,3), 2,1,0, 3,2,1, 4,3,2),

				"b"=c(rep(5,3), 1,2,3, 0,1,2, 1,2,3),
				"d"=c(rep(5,3), 2,1,2, 1,0,1, 2,1,2),
				"g"=c(rep(5,3), 3,2,1, 2,1,0, 3,2,1),


				"m"=c(rep(5,3), 2,3,4, 1,2,3, 0,1,2),
				"n"=c(rep(5,3), 3,2,3, 2,1,2, 1,0,1),
				"q"=c(rep(5,3), 4,3,2, 3,2,1, 2,1,0),

				row.names=sigma
				)

#CLUSTER SETUP

#Initialize empty data frame of clusters; this vector will later change, which will lead to expansion of the constraint set

clusters_matrix <- matrix( nrow = length(sigma), dimnames = list(sigma) )  #"clusters_matrix" is a matrix which has rows named after foci, and every column in the matrix will be a cluster, which is represented as a probability distribution over those foci

#Set up a vector of possible cluster labels; which are combinations of two capital letters

capital_letters <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")

cluster_labels <- do.call( paste, c(expand.grid(capital_letters, capital_letters), sep=""))

cluster_translation <- data.frame( Cluster.name = c(), Reg.expresion = c() ) #Data frame which translates cluster labels to regular expressions


#CONSTRAINT DEFINITIONS

cons.matrix <- NULL #This will be filled in with constraints as they are induced

weights <- c() #This will contain the weights of these constraints


valence <- c(1,-1) #Allow both negative and positive constraints

	#valence <- c(-1) #Allow only negative constraints

#Define constraints

#first, a data frame of all two-symbol constraints; the last position is filled out as "0"

two.symbol <- expand.grid(Valence=valence, First=c(sigma,"^"), Second=c(sigma,"$"), Third="")

two.symbol <- subset(two.symbol, !(two.symbol$First=="^" & two.symbol$Second=="$") ) #Get rid of the constraints "*^$" and "+^$"

#now, a data frame of all three-symbol constraints

three.symbol <- expand.grid(Valence=valence, First=c(sigma,"^"), Second=sigma, Third=c(sigma,"$"))

#finally, vertically concatenate these two data frames using "rbind"

constraints = rbind( two.symbol, three.symbol )

constraints = transform(constraints, Valence = as.numeric(Valence), First = as.character(First), Second = as.character(Second), Third = as.character(Third)) #Convert all columns into columns of characters instead of factor values.

#MAXIMUM GAIN COMPUTATION

#I stole this function from Robert's code so I could add a gradient as an attribute to the gain function

grad <- function(x) attr(x, "gradient")

"grad<-" <- function(x,value) {
  attr(x, "gradient") <- value
  x
}

#The following function returns the maximum gain value 

max.gain <- function(obs.prob, pred.prob, rep_space, new_constr) { #p is the vector of observed probabilities, q is the vector of probabilities predicted by the current grammar, and "new_constr" is the new constraint whose maximal gain value will be computed

	c.star <- violate(new_constr, rep_space) #compute the violation/reward vector of "new_constr" (the new constraint who gain value is to be computed)

	gain <- function(w.star) { #a function of the weight of the new constraint

	gain.value <- w.star * sum(obs.prob*c.star) - log( sum( exp(w.star*c.star) * pred.prob) ) #minus signs are removed from Wilson (2010) definition, and Della Pietra's definition is followed instead

	gain.value <- -gain.value #output the negative of the function, because optim minimizes rather than maximizes functions

	grad(gain.value) <- - ( sum(obs.prob * c.star) - (sum( c.star*exp(w.star*c.star)*pred.prob ) / sum( exp(w.star*c.star)*pred.prob ) )) #Della Pietra's first order derivative

	gain.value #Output the function with the gradient as an attribute to it

	}

	opt.gain <- optim( par = 0, fn = gain, lower = 0, method="L-BFGS-B" ) #minimize the "gain" function with the initial value of w.star being 0

	-opt.gain$value #Output the negative of the value of "gain" associated with its optimal weight (because "gain" was made negative to allow it to be minimized instead of maximized).

}


#Now, a function which determines the violations for a constraint

#Constraints are a series of cells in a data frame!



violate <- function(constr, rep.space) {

	constr.rgxp <- paste(constr[2:4], collapse="") #Bring the characters of the constraint together as a regular expression

	matches.vector <- sapply( rep.space, function(cand) { #omega is the full set of candidates in the representational space

		match <- gregexpr(constr.rgxp, cand)[[1]]

		length( match[match!=-1] )  #Find the number of matches of "const.rgxp" in the candidate; the "match != -1" condition is there because gregexpr returns -1 when there is no match

		}) #This generates a vector of the number of matches for every candidate in omega

	matches.vector*as.numeric(constr[1]) #multiply the vector of matches by 1 or -1, depending on whether the constraint is positive or negative

}



#MAXIMUM ENTROPY CALCULATION

#Compute observed probabilities

p <- sapply( omega, function(cand){ 

				length(which(observed==cand))/length(observed)

			} ) # this is the observed probability (the number of times an element of omega shows up in the observed data divided by the total number of observations)

q <- rep(1/length(omega), length(omega)) # this will be the predicted probability of every form; initially, it will be set to equal probability to every form

prob.frame <- data.frame(omega=omega, p=p, q=q, stringsAsFactors=F) #this is the data frame which will hold p and q for 

#cons.matrix will be a matrix in which each column is the violation vector of a constraint; columns will be indexed with the row number corresponding to that constraint in the "constraints" data frame

#"weights" will be a vector of the weights associated with every constraint in cons.matrix; this vector is provided by the output of the "optim" function 

#(Both are defined in neighbor_search_6.R)

#when new constraints are added, as zeroes are added to "weights" as there are new constraints
#the result is fed to "optim" as initial weights for optimization

#Cover function for optimizing the objective function; biases can be turned on or off.
solve <- function(violation.matrix, wts, lower, var=10000) {

	w <- wts

	obj <- objective.function(violation.matrix)

	obj <- l2.prior(obj,var,mean=0)

	opt <- optim(w, obj, lower=lower, method="L-BFGS-B")

	opt

}


objective.function <- function(violations) { #plug the cons.matrix matrix into this

	objective <- function(w) {

	    	scores <- exp(c(violations %*% w))		#compute maxent scores

          	Z <- sum(scores)

		num <- sum(scores[which(prob.frame$p > 0)])
          
		res <- -(log(num) - log(Z)) #this is the max-ent value for the current weights

	      denom.p <- scores / sum(scores)
      	numer.p <- scores[which(prob.frame$p>0)] / sum(scores[which(prob.frame$p>0)])
        	denom <- colSums(denom.p * violations)
        	numer <- colSums(numer.p * violations[which(prob.frame$p>0), , drop=FALSE])

        	grad(res) <- -(numer - denom) #this is the multi-dimensional gradient for the max-ent function given the current weights
		
		res
	}

	objective

}

# Puts an L2 (log-space Gaussian) prior on an optimization
#Taken from Robert's script
l2.prior <- function(fun, var=1, mean=0, ...) {
  fun <- match.fun(fun)
  function(pv) {

    res <- fun(pv) + (sum(pv - mean)^2) / (2 * var)
    
    if ( !is.null(grad(res)) ) grad(res) <- grad(res) + ((pv - mean) / var)
    
    res
  }
}






#CONSTRAINT SPACE SEARCH

#First, define 3 functions which make neighbor search possible: one to find phonetically close segments, one to create new constraints, and one to find the immediate neighbors of a constraint.

	#function which finds segments which are phonetically close to a certain given segment

phon_dist_search <- function(source_segm, dist_frame) {

	if (source_segm %in% sigma & nchar(source_segm)==1) { #only compute phonetically close values for slots that are not word boundary marks or feature values

		distances <- subset(dist_frame, select = source_segm) 	#find all the distances from a certain segment by selecting the column that is named after that segment
		prox_segments <- rownames(subset(distances,(distances <= phon_dist_threshold) & distances > 0)) #find the names of the rows in which the value in the source segment column is smaller than the maximal proximity threshold

		prox_segments	#output that vector of segments

	}
		
}


	#this function spits out the immediate neighbors of a constraint

immediate_neighbors <- function(constraint_frame, input_row_number, dist_frame) {

	constraint_chars <- as.character(constraint_frame[input_row_number, 1:4])
		
	proximity_list <- lapply(constraint_chars[-1], phon_dist_search, dist_frame)	#create a list of phonetically proximal segments for every character in the constraint


	neighbor_vector <- c()	#initialize the vector if immediate neighbors that is eventually to become the output

	for (i in 1:length(proximity_list)) {	#for every part of the proximity list (corresponding to positions in a constraint), create a vector of constraints derived by replacing the corresponding position in the prototype constraint by each member of the vector of segments created by the phonetic proximity function

		if (length(proximity_list[[i]]) > 0) {	#only if that part of the list is not empty
			for (j in 1:length(proximity_list[[i]]) ) {
			
				new_constraint <- constraint_chars
				new_constraint[i+1] <- proximity_list[[i]][j]	#substitute the character whose turn it is to be taken from the proximity list, into the right position of the original constraint
				output_row_number <- rownames( constraint_frame[ which(      #find the row of the constraints data frame that matches this constraint
													constraint_frame[,1] == new_constraint[1] &    #there should be a more compact way of doing this than just matching every separate column to each separate element of the vector representing the new constraint
													constraint_frame[,2] == new_constraint[2] &
													constraint_frame[,3] == new_constraint[3] &
													constraint_frame[,4] == new_constraint[4]
													),] )
				neighbor_vector <- c(neighbor_vector, output_row_number)
			}

		}

	}

	neighbor_vector

}

#Now comes the part in which these procedures are actually called.


#initialize empty vector of "trash" (=row numbers corresponding to constraints that have been tried and shouldn't be tried again)

trash <- c()

#initialize empty vector of clusters and intersections, to be encoded as regular expressions

total_regexp <- c()

################################ START CYCLE ###################################

#From this point on, the code will be repeated a given number of times

for (iteration in 1:cycles) {   #Repeat all of the following as many times as specified in "cycles"

writeLines(paste("\nCycle number ",iteration,".","\n", sep=""))

#RANDOM CONSTRAINT SELECTION

#Create a sample of the representational space to compute gain values over, to save time
#Currently, the sample is 1 of the rep. space

sample.size <- 1   #size of the sample as a proportion of the total size of omega

sample.indices <- sample( 1:length(omega), round( sample.size * length(omega)), replace=F ) #take a sample of the positions of elements in omega

sample.p <- p[sample.indices]/sum(p[sample.indices])

sample.q <- q[sample.indices]/sum(q[sample.indices])

sample.omega <- omega[sample.indices]

#This repeat loop will ensure that the algorithm keeps going until it finds a row number in the constraint frame which corresponds to a high enough gain value

repeat {

	selected_row <- sample(1:length(constraints$Valence),1) #choose a random row in the constraint frame 

	if ( !(selected_row %in% trash) ) { #If the selected row isn't in the constraints already considered, then do the following; otherwise, go back to repeating

		if ( max.gain(sample.p, sample.q, sample.omega, constraints[selected_row,]) >= accuracy_threshold) {							#if the gain value isn't good enough, then keep trying until you hit a constraint with a sufficiently high gain value
			prototype <- selected_row #if the constraint is accurate enough, declare its original role number to be the prototype for neighbor search
			break    #the while statement is in principle an infinite loop; it is halted whenever a fitting prototype is found
		} else {
			trash <- c(trash, selected_row) #if the constraint isn't accurate enough, add its row number in the original "constraints" frame to "trash"
		}

	}

}

writeLines("Constraint selected at the current cycle:\n")
print(constraints[prototype,])
writeLines("\n")

	#The prototype isn't on the "trash" vector - it hasn't been "used up" yet: it will be used for neighbor search.




#NEIGHBOR SEARCH

## Now, submit this prototype to the neighbor search procedure.

#Initialize a "stack", "trash" and "output_vector" vector.

#trash <- c(0) #"trash" is a vector of all constraint row numbers whose neighbors have been explored

stack <- prototype #"stack" is a vector of constraint row numbers whose neighbors should be explored

trash <- c(trash,prototype) #the prototype should be in "trash" so it doesn't get considered as a neighbor

output_vector <- prototype #"output_vector" is the vector of all constraints found at this iteration of the search algorithm

while (length(stack) > 0 ) {

	current_constraint <- stack[1] #take the first constraint from the stack of constraints to consider

	stack <- stack[-1] #remove the constraint you're considering now from the stack

	neighbors <- immediate_neighbors(constraints, current_constraint, phonetic_distance)	#find the set of immediate neighbors for the prototype
	
	neighbor_frame <- constraints[neighbors,] #look up the shape of these neighboring constraints

	neighbor_frame$Gain <- apply( neighbor_frame, 1, function(x) { max.gain(sample.p, sample.q, sample.omega, x) } ) #compute the gain value of these constraints

	neighbor_frame_pruned <- subset(neighbor_frame, (neighbor_frame$Gain >= accuracy_threshold & !(rownames(neighbor_frame) %in% trash) ))	#select only those constraints that have a sufficient accuracy and which are not among the used constraints
	
	neighbors_pruned <- as.numeric(rownames(neighbor_frame_pruned))			#get the row names of this frame

	output_vector <- c(output_vector,neighbors_pruned) #add the constraint 

	stack <- c(stack, neighbors_pruned) #add the constraints you just found to the stack of constraints whose neighbors should be explored

	trash <- c(trash, neighbors_pruned) #add the constraints you just found to the trash so they don't get found as neighbors again

}

writeLines("Relevant neighboring constraints:\n")
print(constraints[output_vector,])
writeLines("\n")

# put everything in a focus-by-context table

#First, determine all the contexts, and put these in as the first 4 columns of a data frame.

to_be_tabulated <- constraints[output_vector,1:4] #this is the frame of all constraints that were selected by the neighbor search algorithm

# make all possible replacements of non-zero elements with a placeholder "_"

contexts_frame <- data.frame("Valence" = c(), "First" = c(), "Second" = c(), "Third" = c())

for (cons in 1:length(to_be_tabulated[,1])) {	#for every constraint in the "to be tabulated" data frame

	if (to_be_tabulated[cons ,4] == "") {longest <- 3} else {longest <- 4} #if the constraint ends in "0", then only consider the first two positions for replacement, otherwise consider all three trigram positions for replacement

	for (j in 2:longest) { #for every position in the constraint

		context <- to_be_tabulated[cons,]  #take out the row of the contexts_frame that's to be changed now
		if ((context[ , j] != "^") & (context[ , j] != "$") ) { #if the current symbol of the constraint is not a word boundary
			context[ , j] <- "_"		#replace the constraint with "_" in the appropriate position
			contexts_frame <- rbind(contexts_frame, context) #add it to the "contexts_frame" data frame
		}

	}

}

#Get rid of duplicates, so that the list of contexts is as short as possible

contexts_frame <- subset(contexts_frame, !duplicated(contexts_frame)) #this gets rid of all the rows which are the exact same as some row above them

#Now, fill out the focus-by-context data frame by substituting every possible segment into the blanks and computing the corresponding gain value


for (segment in sigma) {     #for every possible segment, put it in the gap left in the context and look up the corresponding gain value

	column_number <- length(contexts_frame) + 1 #this is the number of the column that needs to be added to the contexts frame right now

	focus.substitution.frame <- contexts_frame[,1:4]

	focus.substitution.frame[,2:4] <- apply(focus.substitution.frame[,2:4], 2, function(x){sub("_",segment,x)}) #substitute the current segment value instead of "_" in each context of the contexts_frame

	new_column <- apply( focus.substitution.frame, 1, function(x){ max.gain(sample.p, sample.q, sample.omega, x) } ) #compute the gain value for each constraint made by filling in every context by the focus under scrutiny now

	contexts_frame[,column_number] <- new_column
	names(contexts_frame)[column_number] <- segment

}

#print(contexts_frame)








#CLUSTERING


#The "find_cluster" function applies to a row of "contexts_frame" to find the content of the cluster with the higher mean

find_cluster <- function(rownumber) {

	to_be_clustered <- contexts_frame[rownumber,-(1:4)] #this is the row of the contexts frame over which the one-dimensional clustering analysis will be performed - and the non-numerical values at the left hand side of the frame are excluded

	noise_vector <- sample(c(0 , 0.0000001), length(to_be_clustered), replace = T) #make a random vector of the same length as "to be clustered", with either 0 or 0.0000001 as values

	#This is done so that the clustering algorithm will actually find two clusters in a case like "0 0 0 0 0 0 0 0 0 1 1 1".
	
	to_be_clustered <- to_be_clustered + noise_vector #the noise vector is added to the "to_be_clustered" values (added and not subtracted, so that negative values are not created)

#obtain the parameters for the two Gaussians that were fit, so that these can be submitted to the expectation maximization algorithm

	if ( !is.na(mclustBIC(to_be_clustered,G=1:2,modelNames="E")[2,]) ) { #If it is at all possible to fit a 2-Gaussian model to the data from a given row, then go ahead and find the following values:

		row_Mclust <- Mclust( 	data = to_be_clustered, #take all the columns except the first four
						G = 2, #look for 2 components
						modelNames = "E" ) #use a model in which the variance of both components has to be the same

		if ( length(which(row_Mclust$classification==2)) >= 2 ) {	#if, for the cluster with the higher mean, there is are at least two 

			row_parameters <- row_Mclust$parameters

			#Given these parameter values, perform expectation minimalization 

			row_em <- em(	modelName = "E",
						data = to_be_clustered,
						parameters = row_parameters	)

			classification_matrix <- row_em$z 	#This is a matrix with conditional probabilities that a certain constraint will fall into one of the two clusters

			cluster_content <- classification_matrix[,2] #The content of a cluster is the probability distribution over foci that's associated with the Gaussian mixture component with the higher mean

			return(cluster_content)

		} #If the "cluster" contained just one element, then don't return anything

	} #Endif with the condition that it's possible to fit a two-component model

} #endfunction

#Now, call this function for every row number in the "contexts_frame" frame

all_clusters <- lapply(1:length(contexts_frame$First), find_cluster)

#And after this, get rid of all the members of this list which are NULL.

to_delete <- c()

for (member in 1:length(all_clusters) ) { if (is.null(all_clusters[[member]])) {to_delete <- c(to_delete, member)} } #Find all the members of the all_clusters 

if (length(to_delete) > 0) { all_clusters <- all_clusters[-to_delete] }

#print(all_clusters)



#Procedure for integrating these clusters with previously established clusters (unfinished)

#Then take these clusters, see if they've been established as clusters before, and if not, and if they're "significantly" different from other established clusters, then give them a unique label and add them to the list of clusters

if (length(all_clusters) > 0 & all( is.na(clusters_matrix[,1])) ) {   #If no clusters have been added to the clusters matrix yet, and there are clusters to be added.

	clusters_matrix[,1] <- all_clusters[[1]]	

	label <- sample(cluster_labels, 1)

	colnames(clusters_matrix)[1] <- label

	all_clusters <- all_clusters[-1]

	cluster_labels <- cluster_labels[-which(cluster_labels==label)]

}

#now start comparing every cluster in "all clusters" to existing clusters in "clusters_matrix" and add them to the matrix whenever they are sufficiently different; if they are not sufficiently different, merge the current cluster with the cluster in "cluster_matrix" that it is most similar to

cluster_similarity <- function(clust1, clust2) {  #The distance metric between two clusters (vectors of probabilities assigned to members of sigma) suggested in the meeting on April 3rd

	sum( clust1/sqrt(sum(clust1^2)) * clust2/sqrt(sum(clust2^2)) )

}

similarity_threshold <- .9

#Now, integrate all clusters into the clusters_matrix

while (length(all_clusters) > 0) {


		similarities <- apply( clusters_matrix, 2, cluster_similarity, all_clusters[[1]]) #compute the K-L divergence of the various clusters already present in the cluster matrix given the current cluster you're trying to add

		maximal_similarity <- max(similarities) #Find the greatest similarity value between the cluster to be added and any other given cluster already in the clusters matrix

		clusters_matrix <- cbind(clusters_matrix, all_clusters[[1]])	#add the cluster in question as a new column to the clusters matrix

		if (maximal_similarity < similarity_threshold) {	#if the maximal similarity of the new cluster to any other cluster is not large enough
	
			label <- sample(cluster_labels, 1)	#choose a random, unique label for it

			colnames(clusters_matrix)[ncol(clusters_matrix)] <- label #rename the new column to the label just found

			all_clusters <- all_clusters[-1] #remove the first cluster on the "all_clusters" list

			cluster_labels <- cluster_labels[-which(cluster_labels==label)]	#remove the label just assigned from the range of cluster labels which may be assigned
	
		} else { #if the maximal similarity is at or above the threshold

			cluster_to_be_modified <- which(similarities==maximal_similarity) #find the position of the maximally similar cluster in the cluster matrix
			#If there's more than one column in clusters_matrix that the current cluster is maximally similar to, just pick a value; even when the cluster at hand is equally similar to two clusters with different labels.
			pick_one <- sample(length(cluster_to_be_modified),1)

			cluster_to_be_modified <- cluster_to_be_modified[pick_one] #pick one random member of "cluster_to_be_modified"

			old_label <- names(cluster_to_be_modified) #The vector "cluster_to_be_modified" comes with names for each value that correspond to the name of the cluster the values stand for; so, "old_label" is the cluster label associated with the cluster that the currently considered cluster is about to be merged with.

			colnames(clusters_matrix)[ncol(clusters_matrix)] <- old_label

			all_clusters <- all_clusters[-1]	

		}

}

#print(clusters_matrix)

#UPDATE GRAMMAR
#Now, add the constraints in the neighborhood set to the grammar!

num.of.new.cons <- length(output_vector) #This is how many new constraints are about to be added to the grammar

while (length(output_vector) > 0 ) {

	if (!exists("cons.matrix")) { #if cons.matrix has not yet been created

		cons.matrix <- matrix(
						data = violate(constraints[output_vector[1],], rep.space=omega), #take the first constraint in the output vector, and enter the violations of that constraint as the first column in the constraint
						ncol = 1,
						dimnames = list(omega, output_vector[1]) #the rows should be named after the members of omega; the column is named after the number of the constraint in "constraints"
				    		)

	} else {

		cons.matrix <- cbind(cons.matrix, violate(constraints[output_vector[1],], rep.space=omega) ) #add the vector of violations of the constraint under consideration (obtained through looking it up in "constraints") to the constraint matrix

		colnames(cons.matrix)[ncol(cons.matrix)] <- output_vector[1]

	}
	
	output_vector <- output_vector[-1] #remove the first member of "output_vector", which had just been added to the grammar

}

#Optimize the weights of these constraints


#If the variable "weights" already exists, then simply append as much zeroes as there are new constraints to "weights"; so that the initial weights given to the optimization algorithm are the optimal weights for old constraints and zeroes for new constraints
#Otherwise, create "weights" and put as many zeroes in it as there are new constraints.

if (exists("weights")) {weights <- c(weights, rep(0, num.of.new.cons) ) } else { weights <- rep(0, num.of.new.cons) }

#Give lower bounds to constraint weights to the optimization algorithm
lower <- rep(0, length(weights))

optimized <- solve(cons.matrix, weights, lower=lower, var=1000) #Optimize the weights of the constraints currently in the grammar

#print(optimized)

weights <- optimized$par  #read off the optimized weights

q <- exp(c(cons.matrix %*% weights)) / sum( exp(c(cons.matrix %*% weights)) ) #compute the probabilities of the things in q on the basis of the weights just computed

prob.frame$q <- q #Update the predicted probabilities in the master data frame.

writeLines("Constraint weights assigned at this cycle:\n")
print(cbind(constraints[colnames(cons.matrix),], weights))
writeLines("\n")



#UPDATE CLUSTERS

#When the distributions have been updated, they are projected into statements of discrete sets of segments which are part of every label

#This should be done as follows: all columns of "clusters_matrix" that have a certain label are collected, and the average of probabilities is taken; given this result, all segments which have a probability above a certain cutoff value (for instance, .9) are assigned to that cluster; the content of a cluster may change over the course of the simulation

#The set of all possible constraints is then updated accordingly

cluster_names <- colnames(clusters_matrix)[!(duplicated(colnames(clusters_matrix)))] #Find all unique cluster labels in clusters_matrix

#Write a procedure which will average the values of all columns associated with a particular label, and project from the resulting probabilities a set of segments



new_cluster <- c() #This is to be the vector of clusters that were just added at this cycle of the simulation


for (name in cluster_names) {

	homonym_matrix <- clusters_matrix[,which(colnames(clusters_matrix)==name)] #find all the columns in clusters_matrix which correspond to the same label

	if (length(which(colnames(clusters_matrix)==name)) > 1) {
	
		homonym_average <- rowSums(homonym_matrix)/ncol(homonym_matrix) #average over all the columns which correspond to the label in question

	} else {

		homonym_average <- homonym_matrix #if the homonym matrix only has one column, then the average of the matrix is the same as the matrix itself

	}

	name_content <- names(homonym_average[which(homonym_average > 0.6)]) #The segments that belong to a label are those which belong to that label with more than 0.6 probability.

	reg.expression <- paste(name_content,sep="", collapse="") #Put all the segments that belong to the cluster together

	reg.expression <- paste("[",reg.expression,"]",sep="") #Put brackets around them so they form an actually usable regular expression

	if ( length(which(cluster_translation$Cluster.name==name)) > 0 ) { #if the cluster label is already in cluster_translation

		cluster_translation[which(cluster_translation$Cluster.name==name), 2] <- reg.expression #then just overwrite the value for the reg expression associated with that label

	} else {

		cluster_translation <- rbind(cluster_translation, data.frame( Cluster.name = name, Reg.expression = reg.expression ) )

		new_cluster <- c(new_cluster, name) #only clusters that have just been added are part of "new_cluster"

	}

}



#COMPUTE INTERSECTIONS OF FEATURE CLASSES

if ( length(new_cluster) > 0 ) { #If there are any new clusters

	#For every new cluster, find its intersections with all other clusters in cluster_translation

	new_cluster_regexp <- subset(cluster_translation, cluster_translation$Cluster.name %in% new_cluster)$Reg.expression #Find all the regexpressions corresponding to the new clusters

	cluster_regexp <- cluster_translation$Reg.expression #all reg expressions generated up until now

	two_regexp_combinations <- expand.grid(cluster_regexp, new_cluster_regexp, "", stringsAsFactors = F) #this creates all the combinations of 2 regexpressions, whose intersection is to be computed

	three_regexp_combinations <- expand.grid(cluster_regexp, cluster_regexp, new_cluster_regexp, stringsAsFactors = F) #this creates all the combinations of regexpressions, whose intersection is to be computed

	regexp_combinations <- rbind(two_regexp_combinations, three_regexp_combinations)	#These are all combinations of new regexpressions with new & old regexpressions that have 2 or 3 members

	#Procedure for computing the reg
	compute.intersections <- function(rowcontents) {   #Take three reg. expressions

		rowcontents <- as.character(rowcontents) #

		rowcontents.2 <- lapply(rowcontents, function(x){ 		#Make sure the data frame cell is a character, not a factor level
						reg.char <- unlist(strsplit(x, ""))
						reg.char <- reg.char[-c(1,which(reg.char=="]"))]
						})

		#For each reg.expression, get its characters, and chop off the brackets at the begining and end

		#Intersect the first two reg.expressions.
		intersection1 <- intersect(rowcontents.2[[1]], rowcontents.2[[2]]) 

		#Now, if the third reg. expression is non-empty, then compute the intersection of the first two reg.exp and the third reg. exp; otherwise, do nothing
		if (length(reg3char) > 0) { intersection2 <- intersect(intersection1, rowcontents.2[[3]]) } else { intersection2 <- intersection1 }

		#If the resulting intersection is non-null, paste it together as a regular expression
		if (length(intersection2) > 0) {
			result <- paste(intersection2, collapse="", sep="")
			paste("[", result, "]", collapse="",sep="")
			}

	} #endfunction

	all_intersections <- unlist(apply(regexp_combinations, 1, compute.intersections))

	#Find all unique intersections
	unique_intersections <- unique(all_intersections)

	#Remove all intersections that are identical to earlier created clusters....
	unique_intersections <- unique_intersections[-which(unique_intersections %in% total_regexp)]

	#"New feature combinations" are those that are either newly found clusters, or new intersections of features.
	new_feature_combinations <- c(new_cluster_regexp, unique_intersections)

	#record the clusters and cluster intersections just created in the big book of clusters and intersections
	total_regexp <- c(total_regexp, new_feature_combinations) 

#UPDATE CONSTRAINT SPACE

	#For every member of that set, define all new constraints that contain that member

	#Define all the new two-symbol constraints

	new.two.symbol.1 <- expand.grid(Valence=valence, First=new_feature_combinations, Second=c(sigma,"$"), Third="") #Find all two-symbol constraints such that the first symbol is one of the new clusters/features

	new.two.symbol.2 <- expand.grid(Valence=valence, First=c(sigma,"^"), Second=new_feature_combinations, Third="") #Find all two-symbol constraints

	new.two.symbol.3 <- expand.grid(Valence=valence, First=new_feature_combinations, Second=new_feature_combinations, Third="")


	#Define all the new three-symbol constraints
	#I'm not doing this in one "expand.grid" command because I don't want any of the constraint in the original "constraints" file to be repeated.
	#The larger motivation for this is to keep rownames for constraints constraints


	new.three.symbol.1 <- expand.grid(Valence=valence, First=new_feature_combinations, Second=sigma, Third=c(sigma,"$")) #Only first position 

	new.three.symbol.2 <- expand.grid(Valence=valence, First=c(sigma,"^"), Second=new_feature_combinations, Third=c(sigma,"$"))

	new.three.symbol.3 <- expand.grid(Valence=valence, First=c(sigma,"^"), Second=sigma, Third=new_feature_combinations)

	new.three.symbol.4 <- expand.grid(Valence=valence, First=new_feature_combinations, Second=new_feature_combinations, Third=c(sigma,"$"))

	new.three.symbol.5 <- expand.grid(Valence=valence, First=c(sigma,"^"), Second=new_feature_combinations, Third=new_feature_combinations)

	new.three.symbol.6 <- expand.grid(Valence=valence, First=new_feature_combinations, Second=new_feature_combinations, Third=new_feature_combinations)


	#finally, vertically concatenate these two data frames using "rbind"

	constraints <- rbind(constraints, new.two.symbol.1, new.two.symbol.2, new.two.symbol.3, new.three.symbol.1, new.three.symbol.2, new.three.symbol.3, new.three.symbol.4, new.three.symbol.5, new.three.symbol.6)

	constraints <- transform(constraints, Valence = as.numeric(Valence), First = as.character(First), Second = as.character(Second), Third = as.character(Third)) #Convert all columns into columns of characters instead of factor values.

} else { #endif conditional on there being some new clusters

	new_cluster_regexp <- c() #If there's no clusters to be found at this iteration, output a zero new_cluster_regexp
}

writeLines("Total set of phonological classes induced up until now:\n")
print( cluster_translation$Reg.expression )
writeLines("\n")


#Now, loop this back into the random constraint search procedure.


} #This bracket closes the giant for-loop that this whole part of the code is enclosed in, so that there is a loop from updating the constraint space to random constraint selection.
