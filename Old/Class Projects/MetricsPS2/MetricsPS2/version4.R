#Clinton Tepper
#Problem 4
#NOTE: Starts with the script at the bottom, which calls pCollisions, which in turn calls 
# either branchContrib(...) or contribSimple(...). contribSimple applies to a single collision, while
# branchContrib is slightly more general and can take much longer to calculate.


#Table of probability of a collision given n students
#[1,] 1 0.00000000
#[2,] 2 0.08354969
#[3,] 3 0.23665262
#[4,] 4 0.42789642
#[5,] 5 0.61895982
#[6,] 6 0.77799095
#[7,] 7 0.88915426
#[8,] 8 0.95389151
#[9,] 9 0.98465994
#[10,] 10 0.99617327
#[11,] 11 0.99936375
#[12,] 12 0.99994712

#pCollisions
#gets the probability for a number of collisions greater then or equal to cNum
#Inputs: freqVec (frequency vector),sNum (number of sampls),cNum (min collisions)
#Out: probability
pCollisions = function(sNum, freqVec, cNum) {

    pVec = freqVec / sum(freqVec)
    oNum = length(freqVec) #number of outcomes

    #recursivelly iterates through outcomes
    branchInit = matrix(1:length(freqVec), ncol = 1, nrow = length(freqVec))
    if (cNum > 1)
        rval = sum(apply(branchInit, 1, branchContrib, pVec, sNum, cNum, oNum))
    else #uses the faster method below if possible 
    rval = contribSimple(pVec, sNum, oNum)
    return(rval)

}

#Method 1: branchContrib
#gets the probability for a number of collisions greater then or equal to cNum. Uses permutations.
#Inputs: branch (outcomes), pVec (probabiltiy vector), sNum (number of sampls),cNum (min collisions), oNum (length of p-vector)
#Out: probability for cNum collisions given an input branch of outcomes
branchContrib = function(branch, pVec, sNum, cNum, oNum) {
    branchMat <- matrix(branch, nrow = oNum, ncol = length(branch), byrow = TRUE)
    colMat = matrix(1:oNum, nrow = length(pVec), ncol = length(branch)) == branchMat
    collisions = max(apply(colMat, 1, sum)) #sum the number for each row
    if (collisions > cNum) # check for the number of collisions
    {
        return(prod(c(rep(1, length(branch)) * pVec[branch]))) #return the calculated probability    
    }
    if (length(branch) == sNum) { return(0) }
    #correspons to no collisions found

    #append the branch matrix
    branchMat = cbind(branchMat, 1:oNum)

    #return the sum of the contributions
    return(sum(apply(branchMat, 1, branchContrib, pVec, sNum, cNum, oNum)))

}

#Method 2: Combinatoric method
#Gets the probability for a number of collisions greater then 1
#Inputs: pVec (probabiltiy vector), sNum (number of sampls),oNum (length of p-vector)
#Out: probability
contribSimple = function(pVec, sNum, oNum) {
    if (sNum > oNum)
        return(1)
    combMatrix = combn(1:oNum, sNum) # get all unique combinations
    combMatrix = matrix(pVec[combMatrix], nrow = sNum) #lookup the nessecary probabilities
    #multiply the probabilities to get the permutation probability, accounts for the factorial, and takes the complement
    return(1 - factorial(sNum) * sum(apply(combMatrix, 2, prod)))
}

#WARNING- The 1st method (hits>1) program traverses (almost) the full sample set
# Therefore, any values greater than 5 (scales with O(12^n) in the bday problem will 
# take substantial amounts of time to run IF the hits parameter is greater then one

################################################R Script- Entry Point #############################
nVec = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
fVec = c(340297, 319235, 356786, 329809, 355437, 358251, 367934, 387798, 374711, 367354, 351832, 356111)
NumStudents = 1:12 #number of draws (ie children)
hits = 1 #parameter for the number of coincidences (ie 2 for three people sharing the same birthday, 1 for the classic case of two people sharing one birthday)
mAns = matrix(sapply(NumStudents, pCollisions, fVec, hits), nrow = length(samples))
cat("Table of probability of a collision given n students\n")
mAns = cbind(NumStudents, mAns)
print(mAns)