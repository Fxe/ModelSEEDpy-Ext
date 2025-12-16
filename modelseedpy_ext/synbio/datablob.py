import hashlib
import os
from io import BytesIO
from minio import Minio


class DataBlob:

    def __init__(self, data=None, md5: str = None):
        self.data = data
        self.md5 = md5

        # If data is provided but md5 is missing, process the data
        if data is not None and md5 is None:
            self._process_data(data)

    @property
    def size(self):
        return len(self.data) if self.data else 0

    @staticmethod
    def from_file(file_path):
        """
        Read a file and create a data blob with binary data, MD5 hash, and original size.

        Args:
            file_path (str): Path to the file to read

        Returns:
            DataBlob: Returns a new DataBlob instance
        """
        try:
            # Check if file exists
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")

            # Read the entire file as binary data
            with open(file_path, 'rb') as file:
                data = file.read()

            # Calculate MD5 hash
            md5_hash = hashlib.md5()
            md5_hash.update(data)
            md5 = md5_hash.hexdigest()

            # Create and return DataBlob with data and pre-calculated MD5
            return DataBlob(data=data, md5=md5)

        except Exception as e:
            raise Exception(f"Error reading file {file_path}: {str(e)}")

    def _process_data(self, data):
        """Process data that was provided directly (not from file)"""
        if isinstance(data, str):
            # Convert string to bytes
            self.data = data.encode('utf-8')
        elif isinstance(data, bytes):
            self.data = data
        else:
            raise ValueError("Data must be either string or bytes")

        if self.md5 is None:
            self._calculate_md5()

    def _calculate_md5(self):
        """Calculate MD5 hash of the data"""
        if self.data is not None:
            md5_hash = hashlib.md5()
            md5_hash.update(self.data)
            self.md5 = md5_hash.hexdigest()

    def get_info(self):
        """Get information about the data blob"""
        return {
            'size': self.size,
            'md5': self.md5,
            'data_length': len(self.data) if self.data else 0
        }

    def save_to_file(self, output_path, mkdir: bool = True):
        """Save the data blob to a file"""
        if self.data is None:
            raise ValueError("No data to save")

        if mkdir:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)

        with open(output_path, 'wb') as file:
            file.write(self.data)

    def save_to_minio(self, minio_client: Minio, bucket_name, object_name):
        """
        Save the data blob to MinIO object storage.

        Args:
            minio_client: MinIO client instance
            bucket_name (str): Name of the bucket
            object_name (str): Name of the object in the bucket

        Returns:
            dict: Upload result with object info
        """
        if self.data is None:
            raise ValueError("No data to save")

        try:
            # Create a BytesIO stream from the binary data
            data_stream = BytesIO(self.data)

            # Upload to MinIO
            result = minio_client.put_object(
                bucket_name=bucket_name,
                object_name=object_name,
                data=data_stream,
                length=self.size,
                metadata={
                    'md5': self.md5,
                    'original_size': str(self.size)
                }
            )

            return {
                'bucket_name': bucket_name,
                'object_name': object_name,
                'etag': result.etag,
                'size': self.size,
                'md5': self.md5,
                'version_id': result.version_id if hasattr(result, 'version_id') else None
            }

        except Exception as e:
            raise Exception(f"Error uploading to MinIO: {str(e)}")
        finally:
            # Clean up the stream
            if 'data_stream' in locals():
                data_stream.close()
