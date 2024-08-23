

class Minion:

    def __init__(self, mio_client, bucket):
        self.cli = mio_client
        self.bucket = bucket

    def put_file(self, filename, target, content_type="text/csv"):
        with open(filename, 'rb') as fh:
            self.cli.put_object(self.bucket, target, fh, -1,
                                part_size=10 * 1024 * 1024, content_type=content_type)
