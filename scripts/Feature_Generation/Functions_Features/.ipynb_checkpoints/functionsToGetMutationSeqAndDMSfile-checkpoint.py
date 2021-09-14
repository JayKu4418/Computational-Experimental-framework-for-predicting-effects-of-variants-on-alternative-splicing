import subprocess
import logging
#Function to get insertion or deletion or substitution
# Get fastq file and dms file for insertion
def getFilesForInsertionMutation(position,mutbase,seq,DMSdata=[]):
    seqToRet = seq[0:position]+mutbase+seq[position:]
    
    pos_hit=0
    dmsvals = []
    # GO through every dms val
    if DMSdata:
        # These are the dms values to write that are exactly the same as in the WT file
        dmsvals_to_write_Same = [i for i in DMSdata if int(i[0])<position+1]
        dmsvals_newToWrite = [[str(position+i),"-999.0"] for i in range(1,len(mutbase)+1)]
        dmsvals_to_write_Diff = [[str(int(i[0])+len(mutbase)),i[1]] for i in DMSdata if int(i[0])>position]
        dmsvals = dmsvals_to_write_Same + dmsvals_newToWrite + dmsvals_to_write_Diff
        
    return([seqToRet,dmsvals])

# Get fastq file and dms file for deletion
def getFilesForDeletionMutation(position,wtbase,seq,DMSdata=[]):
    seqToRet = seq[0:position]+seq[position+len(wtbase):]
    pos_hit=0
    dmsvals = []
    if DMSdata:
        # These are the dms values to write that are exactly the same as in the WT file
        dmsvals_to_write_Same = [i for i in DMSdata if int(i[0])<position+1]
        # Only start writing value after the deleted base or bases
        dmsvals_to_write_Diff = [[str(int(i[0])-len(wtbase)),i[1]] for i in DMSdata if int(i[0])>position+len(wtbase)]
        dmsvals = dmsvals_to_write_Same + dmsvals_to_write_Diff
        
    return([seqToRet,dmsvals])

# Get fastq file and dms file for substitution
def getFilesForSubstitutionMutation(position,mutbase,seq,DMSdata=[]):
    seqToRet = seq[0:position]+mutbase+seq[position+1:]
    dmsvals = []
    if DMSdata:
        dmsvals_to_write_Same = [i for i in DMSdata if int(i[0])<position+1]
        dmsvals_newToWrite = [[str(position+i),"-999.0"] for i in range(1,len(mutbase)+1)]
        dmsvals_to_write_Diff = [i for i in DMSdata if int(i[0])>position+len(mutbase)]
        dmsvals = dmsvals_to_write_Same + dmsvals_newToWrite + dmsvals_to_write_Diff
        
    return([seqToRet,dmsvals])

# THis function will run through Rsample, Stochastic, ct2dot and then calculate the energy of each structure in the ensemble
def generateEnsemble(tmp,mutID,fastafile,seed,dmsfile=""):
    # Run through Rsample or Partition, stochastic and ct2dot
    if dmsfile:
        dmsornot="Rsample"
        #subprocess.check_output(["Rsample",fastafile,dmsfile,tmp+mutID+"_RsampleOutput_NoMaxDist.pfs"])
        # New Rsample with Dave DMS parameters
        subprocess.check_output(["Rsample",fastafile,dmsfile,tmp+mutID+"_RsampleOutput_NoMaxDist.pfs","--DMS","--max","5.0"])
    else:
        dmsornot="Partition"
        subprocess.check_output(["partition-smp",fastafile,tmp+mutID+"_PartitionOutput_NoMaxDist.pfs"])
    #print(dmsornot)
    #print(tmp)
    #print(mutID)
    #print(seed)
    seed=str(seed)
    subprocess.check_output(["stochastic",tmp+mutID+"_"+dmsornot+"Output_NoMaxDist.pfs", tmp+mutID+"_"+dmsornot+"Output_NoMaxDist_Seed"+seed+".ct", "--ensemble","1000","--seed", seed])
    subprocess.check_output(["ct2dot", tmp+mutID+"_"+dmsornot+"Output_NoMaxDist_Seed"+seed+".ct", "ALL", tmp+mutID+"_"+dmsornot+"Output_NoMaxDist_Seed"+seed+".db","--format", "simple"])
    subprocess.check_output(["efn2-smp", tmp+mutID+"_"+dmsornot+"Output_NoMaxDist_Seed"+seed+".ct", tmp+mutID+"_"+dmsornot+"Output_NoMaxDist_Seed"+seed+"_efn2.txt"])     
            
def convertEnsembleToElementRepresentation(tmp,mutID,seed,dmsornot):
    # open up db structures and get rid of first two header lines
    with open(tmp+mutID+"_"+dmsornot+"Output_NoMaxDist_Seed"+seed+".db") as f:
        db_strucs = [line.strip() for line in f][2:]
    # Create an empty list to store element structures
    #allstructures=[]
    # Open up write file to write element structures and convert each db structure to element structure
    with open(tmp+mutID+"_"+dmsornot+"Output_NoMaxDist_Seed"+seed+"_Elements.txt","w") as gw:
        for db in db_strucs:
            output = subprocess.check_output(["rnaConvert.py","-T","element_string",str(db)],stderr=subprocess.STDOUT)
            #print(output.decode('ascii').strip().split("\n")[1])
            element = output.decode('ascii').strip().split("\n")[1]
            gw.write(element)
            gw.write("\n")
            #allstructures.append(element)

# This function will run RNAFold to get the MFE structure and then calculate the energy of the structure 
def generateMFE(tmp,mutID,fastafile,dmsfile=""):
    # Run through Rsample or Partition, stochastic and ct2dot
    if dmsfile:
        dmsornot="DMSdata"
        subprocess.check_output(["Fold",fastafile,tmp+mutID+"_MFEstructure_NoMaxDist_DMSdata.ct","--DMS",dmsfile,"--MFE"])
    else:
        dmsornot="NoDMSdata"
        subprocess.check_output(["Fold",fastafile,tmp+mutID+"_MFEstructure_NoMaxDist_NoDMSdata.ct","--MFE"])
    subprocess.check_output(["efn2", tmp+mutID+"_MFEstructure_NoMaxDist_"+dmsornot+".ct", tmp+mutID+"_MFEstructure_NoMaxDist_"+dmsornot+"_efn2.txt"])
    # Get the value of the MFE out
    with open(tmp+mutID+"_MFEstructure_NoMaxDist_"+dmsornot+"_efn2.txt") as f:
        enfn_vals = [float(line.strip().split(' ')[6]) for line in f]
    return(enfn_vals[0])
            
def obtainSequenceAndDMSdataBasedOnDoubleMutation(mut,boundary,WTseq,wt_dms_data=[]):
    mutID=mut[0]
    # This is a double mutation
    if int(mut[1]) < 0:
        pos1 = boundary+int(mut[1])
    else:
        pos1 = boundary-1+int(mut[1])
    if int(mut[4]) < 0:
        pos2 = boundary+int(mut[4])
    else:
        pos2 = boundary-1+int(mut[4])
    # Also switch up bases based on which mutation comes first in position
    if pos1 > pos2:
        position1 = pos2
        WTbase1= mut[5]
        MUTbase1= mut[6]
        position2 = pos1
        WTbase2= mut[2]
        MUTbase2= mut[3]
    else:
        position1 = pos1
        WTbase1= mut[2]
        MUTbase1= mut[3]
        position2 = pos2
        WTbase2= mut[5]
        MUTbase2= mut[6]
    # If first position is an insertion
    if WTbase1=="NA":
        insertiondata=getFilesForInsertionMutation(position1,MUTbase1,WTseq,wt_dms_data)
        seqAfterFirstPos=insertiondata[0]
        dmsvalsAfterFirstPos=insertiondata[1]
        # If an insertion was put in, that means that sequence is 1 longer than it was previously
        # So the position 2 value is 1 off, and 1 must be added to it
        position2=position2+len(MUTbase1)
    # If first position is a deletion
    elif MUTbase1=="NA":
        deletiondata=getFilesForDeletionMutation(position1,WTbase1,WTseq,wt_dms_data)
        seqAfterFirstPos=deletiondata[0]
        dmsvalsAfterFirstPos=deletiondata[1]
        # If an deletion was put in, that means that sequence is 1 shorter than it was previously
        # So the position 2 value is 1 off, and 1 must be removed from it
        position2=position2-len(WTbase1)
    # If first position is neither insertion or deletion
    else:
        if WTbase1==WTseq[position1:position1+len(WTbase1)]:
            substitutiondata=getFilesForSubstitutionMutation(position1,MUTbase1,WTseq,wt_dms_data)
            seqAfterFirstPos=substitutiondata[0]
            dmsvalsAfterFirstPos=substitutiondata[1]
        else:
            logging.debug("Mut " + mutID + " WT bases are not matching bases in WTseq")
            return
    # If second position is an insertion
    if WTbase2=="NA":
        insertiondata=getFilesForInsertionMutation(position2,MUTbase2,seqAfterFirstPos,dmsvalsAfterFirstPos)
        seqToWrite=insertiondata[0]
        dmsvalsToWrite=insertiondata[1]
    # If second position is an deletion
    elif MUTbase2=="NA":     
        deletiondata=getFilesForDeletionMutation(position2,WTbase2,seqAfterFirstPos,dmsvalsAfterFirstPos)
        seqToWrite=deletiondata[0]
        dmsvalsToWrite=deletiondata[1]
    # If second position is neither a deletion nor insertion
    else:
        if WTbase2==seqAfterFirstPos[position2:position2+len(WTbase2)]:
            substitutiondata=getFilesForSubstitutionMutation(position2,MUTbase2,seqAfterFirstPos,dmsvalsAfterFirstPos)
            seqToWrite=substitutiondata[0]
            dmsvalsToWrite=substitutiondata[1]
        else:
            logging.debug("Mut " + mutID + " WT bases are not matching bases in WTseq")
            return
    return([seqToWrite,dmsvalsToWrite])

def obtainSequenceAndDMSdataBasedOnSingleMutation(mut,boundary,WTseq,wt_dms_data=[]):
    mutID=mut[0]
    if int(mut[1]) < 0:
        position = boundary+int(mut[1])
    else:
        position = boundary-1+int(mut[1])
    WTbase = mut[2]
    MUTbase = mut[3]
    # If insertion
    if WTbase=="NA":
        insertiondata=getFilesForInsertionMutation(position,MUTbase,WTseq,wt_dms_data)
        seqToWrite=insertiondata[0]
        dmsvalsToWrite=insertiondata[1]
    # if single deletion
    elif MUTbase=="NA":
        deletiondata=getFilesForDeletionMutation(position,WTbase,WTseq,wt_dms_data)
        seqToWrite=deletiondata[0]
        dmsvalsToWrite=deletiondata[1]
    # If substituion
    elif WTbase==WTseq[position:position+len(WTbase)]: 
    # Check that WTbase matches with the base at the position of the mutation
        subsitutiondata=getFilesForSubstitutionMutation(position,MUTbase,WTseq,wt_dms_data)
        seqToWrite=subsitutiondata[0]
        dmsvalsToWrite=subsitutiondata[1]
    else:
        logging.debug("Mut " + mutID + " WT base not matching base in WTseq")
        return
        
    return([seqToWrite,dmsvalsToWrite])


# Get new start and stop positions of coordinates given double mutations 
def getNewStartAndStopCoordinatesForDoubleMutations(WTbase1,MUTbase1,WTbase2,MUTbase2,pos1,pos2,START,END):
    
    # Determine what type of mutation for first mutations
    
    firstBaseSubstitute=False
    firstBaseInsertion=False
    firstBaseDeletion=False
    
    if WTbase1!="NA" and MUTbase1!="NA":
        firstBaseSubstitute=True
    elif WTbase1!="NA":
        firstBaseDeletion=True
    else:
        firstBaseInsertion=True
        
    # Determine what type of mutation for second mutations
    
    secondBaseSubstitute=False
    secondBaseInsertion=False
    secondBaseDeletion=False
    
    if WTbase2!="NA" and MUTbase2!="NA":
        secondBaseSubstitute=True
    elif WTbase2!="NA":
        secondBaseDeletion=True
    else:
        secondBaseInsertion=True
    
    # only substitutions or first mutation is after the end site 
    if (firstBaseSubstitute and secondBaseSubstitute) or (pos1 > END):
        # no changes 
        newSTART=START
        newEND=END
    # else if first mutation after the start site but before end site, no effect on the start position
    elif pos1 > START:
        newSTART=START
        # if first mutation is substitution
        if firstBaseSubstitute:
            newEND=END
        # if first mutation is insertion
        elif firstBaseInsertion:
            newEND=END + len(MUTbase1)
        # if first mutation is deletion
        elif firstBaseDeletion:
            newEND=END - len(WTbase1)
        # if second mutation after end position, then only first mutation will affect end site 
        # BUT If both first and second mutation within start and end 
        if pos2 <= END:
            # if second mut is insertion
            if secondBaseInsertion:
                newEND=newEND + len(MUTbase2)
            # if second mut is deletion
            elif secondBaseDeletion:
                newEND=newEND - len(WTbase2)
    # else if first mutation before the start site
    elif pos1 <= START:
        # If first mut substitution
        if firstBaseSubstitute:
            newSTART=START
            newEND=END
        # if first mut insertion
        elif firstBaseInsertion:
            newSTART=START + len(MUTbase1)
            newEND=END + len(MUTbase1)
        # If first mut deletion
        elif firstBaseDeletion:
            newSTART=START - len(WTbase1)
            newEND=END - len(WTbase1)
        
        # second mutation is before START
        if pos2 <= END:
            # if second mut is insertion
            if secondBaseInsertion:
                newEND=newEND + len(MUTbase2)
                if pos2 <= START:
                    newSTART=newSTART + len(MUTbase2)
            # if second mut is deletion
            elif secondBaseDeletion:
                newEND=newEND - len(WTbase2)
                if pos2 <= START:
                    newSTART=newSTART - len(WTbase2)
    
    return([newSTART,newEND])