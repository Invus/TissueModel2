def load_subsystem_file():
    with open("allowed_subsystem_list.txt","r") as file:
        bestand=file.readlines()

    #sanitize
    regel_lijst=[]
    for line in bestand:
        new_line=line.replace("\n",'')
        regel_lijst.append(new_line)
    #print(regel_lijst)
    new_line=None
    bestand=None
    file=None
    return regel_lijst

def subsystems():
    allowed_subsystems=load_subsystem_file()
    #read
    print('loading file')
    with open("C:\\Users\\Jelte\\MEP\models\\recon2.2_rep6_fork1.xml","r") as file:
        bestand=file.readlines()
##    with open("C:\\Users\\Jelte\\MEP\models\\tester.txt","r") as file:
##        bestand=file.readlines()
    print('Done')
    new_model=open("C:\\Users\\Jelte\\MEP\models\\recon2.2_minimal12.xml","w")
    
    #sanitize
    print('sanitize and convert to list')
    regel_lijst=[]
    for line in bestand:
        new_line=line.replace("\n",'')
        regel_lijst.append(new_line)
    print('length list: ', len(regel_lijst))
    print('Done')
    new_line=None
    bestand=None
    file=None
    print('parsing file')
    zone=None
    inReaction=False
    lijst={}
    for nr,line in enumerate(regel_lijst):
        if '<listOfReactions>' in line:
            zone='reactions'
            new_model.write(line+'\n')
        elif '</listOfReactions>' in line:
            zone='post_reactions'
            new_model.write(line+'\n')
        elif zone=='reactions' and '<reaction' in line:
            inReaction=True
            subsystem_correct=None
            found_subsystem=False
            #get Id
            start=line.find('id="')+4
            end=line.find('"',start)
            identifier=line[start:end]
            #print(identifier)
            first_line_reaction=nr

        elif inReaction==True and zone=='reactions' and 'SUBSYSTEM:' in line:
            start=line.find('SUBSYSTEM:')+11
            end=line.find('<',start)
            subsystem=line[start:end]
            found_subsystem=True
            
            
        elif inReaction==True and '</reaction>' in line:
            lijst[identifier]=subsystem
            #print(identifier,': ',subsystem)
            inReaction=False
            last_line_reaction=nr
            if found_subsystem==False:
                subsystem='other'
            if subsystem in allowed_subsystems:
                for wrline in regel_lijst[first_line_reaction:last_line_reaction+1]:
                    new_model.write(wrline+'\n')
            
                                
        elif zone!='reactions':
            new_model.write(line+'\n')
        else:
            pass



    
##    #create overview of exsisting subsystems
##    unieke_lijst=[]
##    for value in lijst.values():
##        if value not in unieke_lijst:
##            unieke_lijst.append(value)
##    subsystem_list=open('subsystem_list.txt','w')
##    unieke_lijst.sort()
##    for item in unieke_lijst:
##        subsystem_list.write(item+'\n')
##    subsystem_list.close()
##    
##    print('Done, script completed')
    return lijst

print('starting script')
subsystems()
print('script completed')
