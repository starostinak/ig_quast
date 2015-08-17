from data_structs import ConnectionType

def write_dot_header(handler):
    handler.write('digraph g {\n')
    handler.write('labeldistance=1\n')

def write_dot_ending(handler):
    handler.write('}')

def get_glued_singleton_name(glued_singleton):
    return 's.' + str(glued_singleton.conn_type) + '.' + str(glued_singleton.rep_id) + \
        '.' + str(glued_singleton.cluster_id)

def write_cluster_in_dot(handler, rep_id, cluster_id, cluster_size):
    colors = ['skyblue', 'coral', 'lavender', 'mintcream', 'orange']
    handler.write('"' + str(rep_id) + '.' + str(cluster_id) + '" [label="' + str(cluster_id) + 
                  '\\nsize ' + str(cluster_size) + '" fillcolor=' + colors[rep_id] + ' style=filled]\n')

def write_edge_in_dot(handler, edge):
    handler.write('"' + str(edge.begin[0]) + '.' + str(edge.begin[1]) + '" -> "' + 
        str(edge.end[0]) + '.' + str(edge.end[1]) + '" [label=' + str(len(edge.reads)) +
        ' style=' + ('solid' if edge.conn_type == ConnectionType.distant else 'bold') +
        ' color=' + ('black' if edge.conn_type == ConnectionType.distant else 'red') +
        ('' if edge.conn_type == ConnectionType.close else ' dir=none') + ']\n')

def write_glued_singleton_cluster_in_dot(handler, glued_singleton):
    colors = ['skyblue', 'coral', 'lavender', 'mintcream', 'orange']
    handler.write('"' + get_glued_singleton_name(glued_singleton) + 
        '" [label="' + str(glued_singleton.clusters_count) + \
        ' x size 1" fillcolor=' + colors[0 if glued_singleton.rep_id else 1] + ' style="filled,dashed"]\n')

def write_edge_to_glued_singleton_in_dot(handler, glued_singleton):
    handler.write('"' + get_glued_singleton_name(glued_singleton) + '" -> "' + 
        str(glued_singleton.rep_id) + '.' + str(glued_singleton.cluster_id) + '" [label=' +
        str(glued_singleton.clusters_count) + 
        ' style=' + ('solid' if glued_singleton.conn_type == ConnectionType.distant else 'bold') + 
        ' color=' + ('black' if glued_singleton.conn_type == ConnectionType.distant else 'red') + 
        ('' if glued_singleton.conn_type == ConnectionType.close else ' dir=none') + ']\n')
