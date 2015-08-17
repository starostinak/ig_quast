class ConnectionType:
    close = 1
    close_each = 2
    distant = 3

class Connection:
    def __init__(self, begin, end, reads, conn_type):
        self.begin = begin
        self.end = end
        self.reads = reads
        self.conn_type = conn_type
