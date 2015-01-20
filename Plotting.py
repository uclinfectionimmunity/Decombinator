## plotting functions for decombinator called
## by decombinator_v2.py when withplots == True

def plot_v_usage( handle, chain, species, savefilename="Vusage"):

    ## PLOTS V GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter
    
    if species == 'human':
        
        tags_va = open("humantags_trav.txt", "rU")
        num_genes_va = 0
        gene_list_va = []
        for line in tags_va:
            num_genes_va += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_va.append(x)
        tags_va.close()
        
        tags_vb = open("humantags_trbv.txt", "rU")
        num_genes_vb = 0
        gene_list_vb = []
        for line in tags_vb:
            num_genes_vb += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vb.append(x)
        tags_vb.close()
        
        tags_vg = open("humantags_trgv.txt", "rU")
        num_genes_vg = 0
        gene_list_vg = []
        for line in tags_vg:
            num_genes_vg += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vg.append(x)
        tags_vg.close()
        
        tags_vd = open("humantags_trdv.txt", "rU")
        num_genes_vd = 0
        gene_list_vd = []
        for line in tags_vd:
            num_genes_vd += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vd.append(x)
        tags_vd.close()
            
    elif species == 'mouse':
        
        tags_va = open("mousetags_trav.txt", "rU")
        num_genes_va = 0
        gene_list_va = []
        for line in tags_va:
            num_genes_va += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_va.append(x)
        tags_va.close()
        
        tags_vb = open("mousetags_trbv.txt", "rU")
        num_genes_vb = 0
        gene_list_vb = []
        for line in tags_vb:
            num_genes_vb += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vb.append(x)
        tags_vb.close()
        
        tags_vg = open("mousetags_trgv.txt", "rU")
        num_genes_vg = 0
        gene_list_vg = []
        for line in tags_vg:
            num_genes_vg += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vg.append(x)
        tags_vg.close()
        
        tags_vd = open("mousetags_trdv.txt", "rU")
        num_genes_vd = 0
        gene_list_vd = []
        for line in tags_vd:
            num_genes_vd += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vd.append(x)
        tags_vd.close()
        
    if chain=="alpha":

        freq_vector_v = [0]*num_genes_va
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1
      
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_v)
        percent_usage_v = [0]*num_genes_va
        for i in range(num_genes_va):
            percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
        v_linked = [0]*len(percent_usage_v)
        for i in range(len(percent_usage_v)):
            v_linked[i] = (gene_list_va[i], percent_usage_v[i])
        sorted_v = sorted(v_linked, key=itemgetter(1))
        v_labels = [0]*len(sorted_v)
        v_percents = [0]*len(sorted_v)
        for j in range(len(sorted_v)):
            v_labels[j] = sorted_v[j][0]
            v_percents[j] = sorted_v[j][1]
        pos_v = np.arange(num_genes_va)+ 1
        plt.figure()
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.yticks( pos_v, v_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)
            
    if chain=="beta":

        freq_vector_v = [0]*num_genes_vb
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1
            
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_v)
        percent_usage_v = [0]*num_genes_vb
        for i in range(num_genes_vb):
            percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
        v_linked = [0]*len(percent_usage_v)
        for i in range(len(percent_usage_v)):
            v_linked[i] = (gene_list_vb[i], percent_usage_v[i])
        sorted_v = sorted(v_linked, key=itemgetter(1))
        v_labels = [0]*len(sorted_v)
        v_percents = [0]*len(sorted_v)
        for j in range(len(sorted_v)):
            v_labels[j] = sorted_v[j][0]
            v_percents[j] = sorted_v[j][1]
        pos_v = np.arange(num_genes_vb)+ 1
        plt.figure()
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.yticks( pos_v, v_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="gamma":

        freq_vector_v = [0]*num_genes_vg
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1
 
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_v)
        percent_usage_v = [0]*num_genes_vg
        for i in range(num_genes_vg):
            percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
        v_linked = [0]*len(percent_usage_v)
        for i in range(len(percent_usage_v)):
            v_linked[i] = (gene_list_vg[i], percent_usage_v[i])
        sorted_v = sorted(v_linked, key=itemgetter(1))
        v_labels = [0]*len(sorted_v)
        v_percents = [0]*len(sorted_v)
        for j in range(len(sorted_v)):
            v_labels[j] = sorted_v[j][0]
            v_percents[j] = sorted_v[j][1]
        pos_v = np.arange(num_genes_vg)+ 1
        plt.figure()
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.yticks( pos_v, v_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="delta":

        freq_vector_v = [0]*num_genes_vd
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1
       
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_v)
        percent_usage_v = [0]*num_genes_vd
        for i in range(num_genes_vd):
            percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
        v_linked = [0]*len(percent_usage_v)
        for i in range(len(percent_usage_v)):
            v_linked[i] = (gene_list_vd[i], percent_usage_v[i])
        sorted_v = sorted(v_linked, key=itemgetter(1))
        v_labels = [0]*len(sorted_v)
        v_percents = [0]*len(sorted_v)
        for j in range(len(sorted_v)):
            v_labels[j] = sorted_v[j][0]
            v_percents[j] = sorted_v[j][1]
        pos_v = np.arange(num_genes_vd)+ 1
        plt.figure()
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.yticks( pos_v, v_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    
def plot_j_usage( handle, chain, species, savefilename="Jusage"):

    ## PLOTS J GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter
    
    if species == 'human':

        tags_ja = open("humantags_traj.txt", "rU")
        num_genes_ja = 0
        gene_list_ja = []
        for line in tags_ja:
            num_genes_ja += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_ja.append(x)
        tags_ja.close()
        
        tags_jb = open("humantags_trbj.txt", "rU")
        num_genes_jb = 0
        gene_list_jb = []
        for line in tags_jb:
            num_genes_jb += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jb.append(x)
        tags_jb.close()
        
        tags_jg = open("humantags_trgj.txt", "rU")
        num_genes_jg = 0
        gene_list_jg = []
        for line in tags_jg:
            num_genes_jg += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jg.append(x)
        tags_jg.close()
        
        tags_jd = open("humantags_trdj.txt", "rU")
        num_genes_jd = 0
        gene_list_jd = []
        for line in tags_jd:
            num_genes_jd += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jd.append(x)
        tags_jd.close()
        
    elif species == 'mouse':

        tags_ja = open("mousetags_traj.txt", "rU")
        num_genes_ja = 0
        gene_list_ja = []
        for line in tags_ja:
            num_genes_ja += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_ja.append(x)
        tags_ja.close()
        
        tags_jb = open("mousetags_trbj.txt", "rU")
        num_genes_jb = 0
        gene_list_jb = []
        for line in tags_jb:
            num_genes_jb += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jb.append(x)
        tags_jb.close()
        
        tags_jg = open("mousetags_trgj.txt", "rU")
        num_genes_jg = 0
        gene_list_jg = []
        for line in tags_jg:
            num_genes_jg += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jg.append(x)
        tags_jg.close()
        
        tags_jd = open("mousetags_trdj.txt", "rU")
        num_genes_jd = 0
        gene_list_jd = []
        for line in tags_jd:
            num_genes_jd += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jd.append(x)
        tags_jd.close()

    if chain=="alpha":
           
        freq_vector_j = [0]*num_genes_ja
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_j)
        percent_usage_j = [0]*num_genes_ja
        for i in range(num_genes_ja):
            percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
        j_linked = [0]*len(percent_usage_j)
        for i in range(len(percent_usage_j)):
            j_linked[i] = (gene_list_ja[i], percent_usage_j[i])
        sorted_j = sorted(j_linked, key=itemgetter(1))
        j_labels = [0]*len(sorted_j)
        j_percents = [0]*len(sorted_j)
        for j in range(len(sorted_j)):
            j_labels[j] = sorted_j[j][0]
            j_percents[j] = sorted_j[j][1]
        pos_j = np.arange(num_genes_ja)+ 1
        plt.figure()
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.yticks( pos_j, j_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)
        
    if chain=="beta":
        
        freq_vector_j = [0]*num_genes_jb
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1

        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_j)
        percent_usage_j = [0]*num_genes_jb
        for i in range(num_genes_jb):
            percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
        j_linked = [0]*len(percent_usage_j)
        for i in range(len(percent_usage_j)):
            j_linked[i] = (gene_list_jb[i], percent_usage_j[i])
        sorted_j = sorted(j_linked, key=itemgetter(1))
        j_labels = [0]*len(sorted_j)
        j_percents = [0]*len(sorted_j)
        for j in range(len(sorted_j)):
            j_labels[j] = sorted_j[j][0]
            j_percents[j] = sorted_j[j][1]
        pos_j = np.arange(num_genes_jb)+ 1
        plt.figure()
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.yticks( pos_j, j_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)
        
    if chain=="gamma":

        freq_vector_j = [0]*num_genes_jg
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_j)
        percent_usage_j = [0]*num_genes_jg
        for i in range(num_genes_jg):
            percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
        j_linked = [0]*len(percent_usage_j)
        for i in range(len(percent_usage_j)):
            j_linked[i] = (gene_list_jg[i], percent_usage_j[i])
        sorted_j = sorted(j_linked, key=itemgetter(1))
        j_labels = [0]*len(sorted_j)
        j_percents = [0]*len(sorted_j)
        for j in range(len(sorted_j)):
            j_labels[j] = sorted_j[j][0]
            j_percents[j] = sorted_j[j][1]
        pos_j = np.arange(num_genes_jg)+ 1
        plt.figure()
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.yticks( pos_j, j_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)
        
    if chain=="delta":
        
        freq_vector_j = [0]*num_genes_jd
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_j)
        percent_usage_j = [0]*num_genes_jd
        for i in range(num_genes_jd):
            percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
        j_linked = [0]*len(percent_usage_j)
        for i in range(len(percent_usage_j)):
            j_linked[i] = (gene_list_jd[i], percent_usage_j[i])
        sorted_j = sorted(j_linked, key=itemgetter(1))
        j_labels = [0]*len(sorted_j)
        j_percents = [0]*len(sorted_j)
        for j in range(len(sorted_j)):
            j_labels[j] = sorted_j[j][0]
            j_percents[j] = sorted_j[j][1]
        pos_j = np.arange(num_genes_jd)+ 1
        plt.figure()
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.yticks( pos_j, j_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)
        
    handle.close()

def plot_del_v( handle, savefilename="Vdels"):

    ## PLOTS V GERMLINE DELETIONS BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_v = [0]*50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_v[int(elements.split(',')[2])] += 1

    total = sum(deletions_v)
    for i in range(len(deletions_v)):
        deletions_v[i] = deletions_v[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_v, width, color='yellow')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of V germline deletions', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.ylim((0,0.2))
    plt.xlim((0,20))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()

def plot_del_j( handle, savefilename="Jdels"):

    ## PLOTS J GERMLINE DELETIONS BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_j = [0]*50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_j[int(elements.split(',')[3])] += 1

    total = sum(deletions_j)
    for i in range(len(deletions_j)):
        deletions_j[i] = deletions_j[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_j, width, color='red')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of J germline deletions', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.ylim((0,0.2))
    plt.xlim((0,20))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()

def plot_vj_joint_dist( handle, species, chain, savefilename="VJusage" ):

    ## PLOTS VJ JOINT GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    
    if species == 'human':

        tags_va = open("humantags_trav.txt", "rU")
        num_genes_va = 0
        gene_list_va = []
        for line in tags_va:
            num_genes_va += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_va.append(x)
        tags_va.close()
        
        tags_vb = open("humantags_trbv.txt", "rU")
        num_genes_vb = 0
        gene_list_vb = []
        for line in tags_vb:
            num_genes_vb += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vb.append(x)
        tags_vb.close()
        
        tags_vg = open("humantags_trgv.txt", "rU")
        num_genes_vg = 0
        gene_list_vg = []
        for line in tags_vg:
            num_genes_vg += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vg.append(x)
        tags_vg.close()
        
        tags_vd = open("humantags_trdv.txt", "rU")
        num_genes_vd = 0
        gene_list_vd = []
        for line in tags_vd:
            num_genes_vd += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vd.append(x)
        tags_vd.close()
            
    elif species == 'mouse':
        
        tags_va = open("mousetags_trav.txt", "rU")
        num_genes_va = 0
        gene_list_va = []
        for line in tags_va:
            num_genes_va += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_va.append(x)
        tags_va.close()
        
        tags_vb = open("mousetags_trbv.txt", "rU")
        num_genes_vb = 0
        gene_list_vb = []
        for line in tags_vb:
            num_genes_vb += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vb.append(x)
        tags_vb.close()
        
        tags_vg = open("mousetags_trgv.txt", "rU")
        num_genes_vg = 0
        gene_list_vg = []
        for line in tags_vg:
            num_genes_vg += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vg.append(x)
        tags_vg.close()
        
        tags_vd = open("mousetags_trdv.txt", "rU")
        num_genes_vd = 0
        gene_list_vd = []
        for line in tags_vd:
            num_genes_vd += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_vd.append(x)
        tags_vd.close()
            
    if species == 'human':

        tags_ja = open("humantags_traj.txt", "rU")
        num_genes_ja = 0
        gene_list_ja = []
        for line in tags_ja:
            num_genes_ja += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_ja.append(x)
        tags_ja.close()
        
        tags_jb = open("humantags_trbj.txt", "rU")
        num_genes_jb = 0
        gene_list_jb = []
        for line in tags_jb:
            num_genes_jb += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jb.append(x)
        tags_jb.close()
        
        tags_jg = open("humantags_trgj.txt", "rU")
        num_genes_jg = 0
        gene_list_jg = []
        for line in tags_jg:
            num_genes_jg += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jg.append(x)
        tags_jg.close()
        
        tags_jd = open("humantags_trdj.txt", "rU")
        num_genes_jd = 0
        gene_list_jd = []
        for line in tags_jd:
            num_genes_jd += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jd.append(x)
        tags_jd.close()
            
    elif species == 'mouse':
        
        tags_ja = open("mousetags_traj.txt", "rU")
        num_genes_ja = 0
        gene_list_ja = []
        for line in tags_ja:
            num_genes_ja += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_ja.append(x)
        tags_ja.close()
        
        tags_jb = open("mousetags_trbj.txt", "rU")
        num_genes_jb = 0
        gene_list_jb = []
        for line in tags_jb:
            num_genes_jb += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jb.append(x)
        tags_jb.close()
        
        tags_jg = open("mousetags_trgj.txt", "rU")
        num_genes_jg = 0
        gene_list_jg = []
        for line in tags_jg:
            num_genes_jg += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jg.append(x)
        tags_jg.close()
        
        tags_jd = open("mousetags_trdj.txt", "rU")
        num_genes_jd = 0
        gene_list_jd = []
        for line in tags_jd:
            num_genes_jd += 1
            x = line.split("|")[1].split("*")[0].split("/")[0]
            gene_list_jd.append(x)
        tags_jd.close()
        
    if chain=="alpha":
        
        joint_distribution = np.zeros((num_genes_va,num_genes_ja))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
     
        pos_v = np.arange(num_genes_va)+ 1
        pos_j = np.arange(num_genes_ja)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_va)
        plt.xticks( pos_ticks_j, gene_list_ja)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, rotation='vertical', fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="beta":
        
        joint_distribution = np.zeros((num_genes_vb,num_genes_jb))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))

        pos_v = np.arange(num_genes_vb)+ 1
        pos_j = np.arange(num_genes_jb)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_vb)
        plt.xticks( pos_ticks_j, gene_list_jb)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="gamma":
        
        joint_distribution = np.zeros((num_genes_vg,num_genes_jg))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))

        pos_v = np.arange(num_genes_vg)+ 1
        pos_j = np.arange(num_genes_jg)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_vg)
        plt.xticks( pos_ticks_j, gene_list_jg)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="delta":
        
        joint_distribution = np.zeros((num_genes_vd,num_genes_jd))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
 
        pos_v = np.arange(num_genes_vd)+ 1
        pos_j = np.arange(num_genes_jd)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_vd)
        plt.xticks( pos_ticks_j, gene_list_jd)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()

def plot_insert_lengths( handle, savefilename="InsertLengths" ):

    ## PLOTS DISTRIBUTION OF INSERT LENGTH BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    maxi = 500
    insert_lengths = [0]*maxi

    for line in handle:
        elements = line.rstrip("\n")

        classifier = elements.split(',')
        if len(classifier) == 8:
            insert_lengths[len(classifier[4].replace(' ',''))] += 1
        else:
            insert_lengths[0] += 1

    total = sum(insert_lengths)
    for i in range(len(insert_lengths)):
        insert_lengths[i] = insert_lengths[i] / float(total)

    ind = np.arange(maxi)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, insert_lengths, width, color='blue')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of nucleotides', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.xlim((0,50))

    plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()