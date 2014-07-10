## plotting functions for decombinator called
## by decombinator_v2.py when withplots == True

def plot_v_usage( handle, chain, savefilename="Vusage", order="frequency"):

    ## PLOTS V GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    if chain=="alpha":
        tags = open("tags_trav.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order=="frequency":        
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_v)
            percent_usage_v = [0]*num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V1-1','V1-2','V10','V12-1','V12-2','V12-3','V13-1','V13-2','V14/D4','V16','V17','V18','V19','V2','V20','V21','V22','V23/D6','V24','V25','V26-1','V26-2','V27','V29/DV5','V3','V30','V34','V35','V36/DV7','V38-1','V38-2/DV8','V39','V4','V40','V41','V5','V6','V7','V8-1','V8-2/8-4','V8-3','V8-6','V9-1','V9-2','DV1','DV2','DV3')
            v_linked = [0]*len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0]*len(sorted_v)
            v_percents = [0]*len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.yticks( pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)
            
        elif order=="chromosome":
            total = sum(freq_vector_v)
            fv = [0]*num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V1-1','V1-2','V10','V12-1','V12-2','V12-3','V13-1','V13-2','V14/D4','V16','V17','V18','V19','V2','V20','V21','V22','V23/D6','V24','V25','V26-1','V26-2','V27','V29/DV5','V3','V30','V34','V35','V36/DV7','V38-1','V38-2/DV8','V39','V4','V40','V41','V5','V6','V7','V8-1','V8-2/8-4','V8-3','V8-6','V9-1','V9-2','DV1','DV2','DV3')
            chromosome_order = [0,1,13,24,32,35,36,37,38,42,2,3,39,40,6,4,7,8,43,5,41,9,10,11,12,14,15,16,17,44,18,19,20,22,23,25,21,26,27,28,29,30,31,33,34,45,46]
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize = 8)
            ax.set_xticks(ind+1*width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="beta":
        tags = open("tags_trbv.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order=="frequency":        
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_v)
            percent_usage_v = [0]*num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V10-1','V10-2','V10-3','V11-1','V11-2','V11-3','V12-3/V12-4','V12-5','V13','V14','V15','V16','V18','V19','V2','V20-1','V24-1','V25-1','V27-1','V28-1','V29-1','V3-1','V30-1','V4-1','V4-2','V4-3','V5-1','V5-4','V5-5','V5-6','V5-8','V6-1','V6-4','V6-5','V6-6','V6-8','V6-9','V7-2','V7-3','V7-4','V7-6','V7-7','V7-8','V7-9','V9')
            v_linked = [0]*len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0]*len(sorted_v)
            v_percents = [0]*len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.yticks( pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)
            
        elif order=="chromosome":
            total = sum(freq_vector_v)
            fv = [0]*num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V10-1','V10-2','V10-3','V11-1','V11-2','V11-3','V12-3/V12-4','V12-5','V13','V14','V15','V16','V18','V19','V2','V20-1','V24-1','V25-1','V27-1','V28-1','V29-1','V3-1','V30-1','V4-1','V4-2','V4-3','V5-1','V5-4','V5-5','V5-6','V5-8','V6-1','V6-4','V6-5','V6-6','V6-8','V6-9','V7-2','V7-3','V7-4','V7-6','V7-7','V7-8','V7-9','V9')
            chromosome_order = [ 14, 21, 23, 26, 31, 24, 25, 37, 32, 38, 44, 0, 3, 1, 4, 33, 39, 27, 34, 28, 40, 29, 35, 41, 36, 42, 30, 43, 8, 2, 5, 6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 22 ]
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+1*width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="gamma":
        tags = open("tags_trgv.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order=="frequency":        
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_v)
            percent_usage_v = [0]*num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V2','V3','V4','V5','V8','V9')
            v_linked = [0]*len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0]*len(sorted_v)
            v_percents = [0]*len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.yticks( pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)
            
        elif order=="chromosome":
            total = sum(freq_vector_v)
            fv = [0]*num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V2','V3','V4','V5','V8','V9')
            chromosome_order = [0,1,2,3,4,5]
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+1*width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="delta":
        tags = open("tags_trdv.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order=="frequency":        
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_v)
            percent_usage_v = [0]*num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V1','V2','V3','V4','V5','V6','V7','V8')
            v_linked = [0]*len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0]*len(sorted_v)
            v_percents = [0]*len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.yticks( pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)

        elif order=="chromosome":
            total = sum(freq_vector_v)
            fv = [0]*num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V1','V2','V3','V4','V5','V6','V7','V8')
            chromosome_order = [ 3, 5, 0, 4, 6, 7, 1, 2 ]
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+1*width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    tags.close()

    
def plot_j_usage( handle, chain="beta", savefilename="Jusage", order="frequency"):

    ## PLOTS J GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    if chain=="alpha":
        
        tags = open("tags_traj.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order=="frequency":
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_j)
            percent_usage_j = [0]*num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J10','J11','J12','J13','J14','J15','J16','J17','J18','J20','J21','J22','J23','J24','J26','J27','J28','J29','J3','J30','J31','J32','J33','J34','J36','J37','J38','J39','J4','J40','J41','J42','J43','J44','J45','J46','J47','J48','J49','J5','J50','J52','J53','J54','J56','J57','J6','J7','J8','J9')
            j_linked = [0]*len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0]*len(sorted_j)
            j_percents = [0]*len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.yticks( pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)

        elif order=="chromosome":
            total = sum(freq_vector_j)
            fj = [0]*num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J10','J11','J12','J13','J14','J15','J16','J17','J18','J20','J21','J22','J23','J24','J26','J27','J28','J29','J3','J30','J31','J32','J33','J34','J36','J37','J38','J39','J4','J40','J41','J42','J43','J44','J45','J46','J47','J48','J49','J5','J50','J52','J53','J54','J56','J57','J6','J7','J8','J9')
            chromosome_order = [45,44,43,42,41,40,38,37,36,35,34,33,32,31,30,29,27,26,25,24,23,22,21,20,19,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,49,48,47,46,39,28,18]
            gene_list_j = [ gene_list_j[i] for i in chromosome_order]
            fj = [fj[i] for i in chromosome_order]
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="beta":
        
        tags = open("tags_trbj.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order=="frequency":
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_j)
            percent_usage_j = [0]*num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')
            j_linked = [0]*len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0]*len(sorted_j)
            j_percents = [0]*len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.yticks( pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)

        elif order=="chromosome":
            total = sum(freq_vector_j)
            fj = [0]*num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="gamma":
         
        tags = open("tags_trgj.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order=="frequency":
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_j)
            percent_usage_j = [0]*num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1','J2','JP','JP1','JP2')
            j_linked = [0]*len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0]*len(sorted_j)
            j_percents = [0]*len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.yticks( pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)

        elif order=="chromosome":
            total = sum(freq_vector_j)
            fj = [0]*num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1','J2','JP','JP1','JP2')
            chromosome_order = [3,2,0,4,1]
            gene_list_j = [ gene_list_j[i] for i in chromosome_order]
            fj = [fj[i] for i in chromosome_order]
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="delta":
         
        tags = open("tags_trdj.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0]*num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order=="frequency":
            plt.rcParams['figure.figsize'] = 10,10
            total = sum(freq_vector_j)
            percent_usage_j = [0]*num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1','J2','J3','J4')
            j_linked = [0]*len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0]*len(sorted_j)
            j_percents = [0]*len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes)+ 1
            plt.figure()
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.yticks( pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename)+'.png', dpi=300)

        elif order=="chromosome":
            total = sum(freq_vector_j)
            fj = [0]*num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1','J2','J3','J4')
            chromosome_order = [0,3,1,2]
            gene_list_j = [ gene_list_j[i] for i in chromosome_order]
            fj = [fj[i] for i in chromosome_order]
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize = 10)
            ax.set_xticks(ind+width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 6)
            plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 6)
            plt.grid(True)

            plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    tags.close()

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

def plot_vj_joint_dist( handle, chain="beta", savefilename="VJusage" ):

    ## PLOTS VJ JOINT GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
        
    if chain=="alpha":
        
        tags_v = open("tags_trav.txt", "rU")
        tags_j = open("tags_traj.txt", "rU")

        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v,num_j))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
        gene_list_v = ('V1-1','V1-2','V10','V12-1','V12-2','V12-3','V13-1','V13-2','V14/D4','V16','V17','V18','V19','V2','V20','V21','V22','V23/D6','V24','V25','V26-1','V26-2','V27','V29/DV5','V3','V30','V34','V35','V36/DV7','V38-1','V38-2/DV8','V39','V4','V40','V41','V5','V6','V7','V8-1','V8-2/8-4','V8-3','V8-6','V9-1','V9-2','DV1','DV2','DV3')
        gene_list_j = ('J10','J11','J12','J13','J14','J15','J16','J17','J18','J20','J21','J22','J23','J24','J26','J27','J28','J29','J3','J30','J31','J32','J33','J34','J36','J37','J38','J39','J4','J40','J41','J42','J43','J44','J45','J46','J47','J48','J49','J5','J50','J52','J53','J54','J56','J57','J6','J7','J8','J9')
            
        pos_v = np.arange(num_v)+ 1
        pos_j = np.arange(num_j)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_v)
        plt.xticks( pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, rotation='vertical', fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="beta":
        tags_v = open("tags_trbv.txt", "rU")
        tags_j = open("tags_trbj.txt", "rU")
        
        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v,num_j))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
        gene_list_v = ('V10-1','V10-2','V10-3','V11-1','V11-2','V11-3','V12-3/V12-4','V12-5','V13','V14','V15','V16','V18','V19','V2','V20-1','V24-1','V25-1','V27-1','V28-1','V29-1','V3-1','V30-1','V4-1','V4-2','V4-3','V5-1','V5-4','V5-5','V5-6','V5-8','V6-1','V6-4','V6-5','V6-6','V6-8','V6-9','V7-2','V7-3','V7-4','V7-6','V7-7','V7-8','V7-9','V9')
        gene_list_j = ('J1-1','J1-2','J1-3','J1-4','J1-5','J1-6','J2-1','J2-2','J2-3','J2-4','J2-5','J2-6','J2-7')

        pos_v = np.arange(num_v)+ 1
        pos_j = np.arange(num_j)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_v)
        plt.xticks( pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="gamma":
        tags_v = open("tags_trgv.txt", "rU")
        tags_j = open("tags_trgj.txt", "rU")
        
        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v,num_j))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
        gene_list_v = ('V2','V3','V4','V5','V8','V9')
        gene_list_j = ('J1','J2','JP','JP1','JP2')
        
        pos_v = np.arange(num_v)+ 1
        pos_j = np.arange(num_j)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_v)
        plt.xticks( pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    if chain=="delta":
        tags_v = open("tags_trdv.txt", "rU")
        tags_j = open("tags_trdj.txt", "rU")
        
        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v,num_j))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v,j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
        gene_list_v = ('V1','V2','V3','V4','V5','V6','V7','V8')
        gene_list_j = ('J1','J2','J3','J4')
        
        pos_v = np.arange(num_v)+ 1
        pos_j = np.arange(num_j)+ 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v-0.5
        pos_ticks_j = pos_j-0.5
        plt.yticks( pos_ticks_v, gene_list_v)
        plt.xticks( pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename)+'.png', dpi=300)

    handle.close()
    tags_v.close()
    tags_j.close()

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