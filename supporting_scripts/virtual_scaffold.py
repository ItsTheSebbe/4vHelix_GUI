import random
import itertools
from openpyxl import Workbook
from pathlib import Path

def random_seq_creator(length):

    bases = ["A", "C", "G", "T"]
    sequence = []

    for _ in range(length):
        base = random.choices(bases, weights=(29, 21, 21, 29), k = 1)
        sequence = sequence + (base)

    print(sequence)
    sequence = ''.join(sequence)
    return sequence

def consecutive_g_count(sequence):

    consecutive_G = 0

    counter = ([[k, len(list(g))] for k, g in itertools.groupby(sequence)])
    #G_counter = [char for char in counter if counter [0]=="G"]
    for char in counter:
        if char[0] == "G":
            G_count = char[1]
            if G_count >= consecutive_G:
                consecutive_G = G_count

    return consecutive_G

def consecutive_c_count(sequence):

    consecutive_C = 0

    counter = ([[k, len(list(g))] for k, g in itertools.groupby(sequence)])
    #C_counter = [char for char in counter if counter [0]=="G"]
    for char in counter:
        if char[0] == "C":
            C_count = char[1]
            if C_count >= consecutive_C:
                consecutive_C = C_count

    return consecutive_C

def slicing_sequence(sequence, cut_template, edge_number):

    sliced_seq = {}

    start = 0
    position = 1

    for i in cut_template:
        chunk = sequence[start:start+i]
        start = start + i
        sliced_seq['edge '+ str(edge_number+1) + ': segment ' + str(position) + ' (length = ' + str(i) + ')'] = chunk
        position = position + 1

    return sliced_seq

def gc_content(sequence, length):

    GC_count = 0

    for base in sequence:
        if base == 'G' or base == 'C':
            GC_count = GC_count + 1

    GC_percentage = GC_count/length * 100

    return GC_percentage

def sequence_creator(length):

    should_restart = True

    while should_restart == True:

        should_restart = False
        sequence = random_seq_creator(length)
        consecutive_G = consecutive_g_count(sequence)
        consecutive_C = consecutive_c_count(sequence)
        gc_percentage = gc_content(sequence, length)

        if consecutive_C > 4 or consecutive_G > 4 or gc_percentage > 44:
            should_restart = True
            continue

    #print ('Finished,', "Consecutive C = ", consecutive_C, "Consecutive G = ", consecutive_G)

    return sequence, gc_percentage

def save_workbook(edge_number, bundles, filename, deletions_list, double_vertices, faces_to_edges_list, to_reinforce, double_edges):

    workbook = Workbook()
    sheet = workbook.active
    row = 1
    i = 3

    to_reinforce[:] = [x for x in to_reinforce if x not in double_edges] #remove the double edges from the list to be reinforced https://stackoverflow.com/questions/1207406/how-to-remove-items-from-a-list-while-iterating

    print(to_reinforce)

    for edge in to_reinforce:
        random.seed(edge)       #added a seed for the virtual scaffold sequence.
                                # It uses the number of the edge so the sequence is different between
                                # different edges, but always the same for the same dge

        if double_vertices == True and edge+1 in to_reinforce:

            sequence, gc_percentage = sequence_creator(bundles['bundle_' + str(edge)][0]['helix_0']['length']) #virtual scaff with autofill
            length = bundles['bundle_' + str(edge)][0]['helix_0']['length'] - deletions_list[edge-1]        #not sure why -1, to look into
            length_for_cut = round((length - 18) / 32)
            cut_template = []
            cut_template.append(18)
            cut_template[1:] = [32 for x in range(length_for_cut-1)]


            cut_template.append(int(((bundles['bundle_' + str(edge)][0]['helix_0']['length']-18)/32-length_for_cut)*32 + 32))

            sliced_sequences = (slicing_sequence(sequence, cut_template, edge-1)) #not sure why -1, to look into

            sheet['A'+str(row)] = 'Sequence total edge '+ str(edge)
            sheet['B'+str(row)] = str(gc_percentage) + '%'
            sheet['C'+str(row)] = sequence
            sheet['D'+str(row+1)] = '=len(C'+str(row)+')'
            sheet['B'+str(row+1)] = 'Number of segments edge ' + str(edge) + ' '
            sheet['C'+str(row+1)] = len(sliced_sequences)

            for key, item in sliced_sequences.items():
                sheet['A' + str(i+1)] = key
                sheet['B' + str(i+1)] = item
                sheet['C' + str(i+1)] = '=len(B' + str(i+1) + (')')
                i = i + 1
            i = i + 4
            row = i-2

        else:

            sequence, gc_percentage = sequence_creator(bundles['bundle_' + str(edge)][0]['helix_3']['length'])
            length = bundles['bundle_' + str(edge)][0]['helix_3']['length'] - deletions_list[edge-1] #not sure why -1, to look into
            length_for_cut = round((length - 18) / 32)

            cut_template = []
            cut_template.append(18)
            cut_template[1:] = [32 for x in range(length_for_cut - 1)]

            cut_template.append(int(((bundles['bundle_' + str(edge)][0]['helix_3'][
                                          'length'] - 18) / 32 - length_for_cut) * 32 + 32))

            sliced_sequences = (slicing_sequence(sequence, cut_template, edge)) #not sure why -1, to look into

            sheet['A' + str(row)] = 'Sequence total edge ' + str(edge)
            sheet['B' + str(row)] = str(gc_percentage) + '%'
            sheet['C' + str(row)] = sequence
            sheet['D' + str(row + 1)] = '=len(C' + str(row) + ')'
            sheet['B' + str(row + 1)] = 'Number of segments edge ' + str(edge) + ' '
            sheet['C' + str(row + 1)] = len(sliced_sequences)

            for key, item in sliced_sequences.items():
                sheet['A' + str(i+1)] = key
                sheet['B' + str(i+1)] = item
                sheet['C' + str(i+1)] = '=len(B' + str(i+1) + (')')
                i = i + 1
            i = i + 4
            row = i - 2

    i = 3
    row_target = 1

    if double_vertices == True :

        virt_scaf_cross_previous = 1


        for edge in range(1, edge_number+1):

            for x in faces_to_edges_list:
                if edge in x:
                    virt_scaf_cross_previous = x[(x.index(edge) - 1)]
                    try:
                        virt_scaf_cross_next = x[x.index(edge) + 1]
                    except IndexError:
                        virt_scaf_cross_next = x[0]

            sheet['F' + str(i+1)] = 'Edge ' + str(edge) + '-' +  str(virt_scaf_cross_previous) + ' connected'
            number_of_segments = sheet['C'+str(i-1)].value

            for row in sheet.iter_rows():   #search for the number of segment and coordinate previous strand (the one to connect)
                for ele in row:
                    if ele.value != None and ('Number of segments edge '+ str(virt_scaf_cross_previous) + ' ' ) in str(ele.value): #loop to find the number of segment target and
                        number_of_segments_target = sheet['C' + str(ele.row)].value
                        coordinate_cell_target =  sheet['B'+ str(ele.row + number_of_segments_target-1)].coordinate
                        sheet['G' + str(i+1)] = '=_xlfn.CONCAT(' + str(sheet['B' + str(i+1)].coordinate) + ',' + str(coordinate_cell_target) + ')'
                        sheet['H' + str(i+2)] = '=LEN(' + 'G' + str(i+1) + ')'

            i = i + 4 + number_of_segments

                    #if ele != None and ('edge ' + str(virt_scaf_cross_previous) + ':segment' + str(number_of_segments_target)) in str(ele):
                     #   print(ele)
    #while sheet["A1"].value is not None:
    #    print("ok")

    # filename = os.path.join(folder, filename)
    file_save = "./" + str(filename[:-16]) + "/" + filename
    workbook.save(file_save)

if __name__ == "__main__":

    length = 220

    edge_length = 220
    edge_number = 2

    filename="test.xlsx"
    save_workbook(edge_number, length, filename, double_edges)

