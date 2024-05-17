# Import necessary libraries
from Crypto.Cipher import AES
from Bio.Seq import Seq
from Crypto.Util.Padding import pad, unpad
from Crypto.Random import get_random_bytes
import ipfshttpclient
import os
import binascii
import timeit

# Define a class for DNA Encryption System
class DNAEncryptionSystem:

    # Constructor to initialize plaintext and encryption key
    def __init__(self, my_plaintext_, encryption_key):
        self.my_plaintext_ = my_plaintext_
        self.encryption_key = encryption_key

    # Convert binary string to DNA sequence
    def convert_dna_format(self, _bin_string):
        map_dna = {'00': 'A', '01': 'T', '10': 'G', '11': 'C'}
        
        # Ensure the binary string has an even length
        if len(_bin_string) % 2 != 0:
            _bin_string += '0'
        
        # Map binary pairs to DNA bases
        get_dna_sequence = ''.join([map_dna[_bin_string[i:i+2]] for i in range(0, len(_bin_string), 2)])
        return get_dna_sequence

    # Apply DNA encryption to the encryption key
    def apply_dna_encryption(self):
        # Convert encryption key to binary
        binary_representation = bin(int.from_bytes(self.encryption_key.encode(), 'big'))[2:]
        # Convert binary to DNA sequence and reverse complement
        get_dna_sequence = self.convert_dna_format(binary_representation)
        return str(Seq(get_dna_sequence).reverse_complement())

    # Perform AES encryption with custom padding
    def perform_aes_encryption(self, key, my_plaintext_):
        _blk_size = 16

        # Convert plaintext to bytes
        plaintext_bytes = my_plaintext_.encode()

        # Apply custom padding
        padded_plaintext = pad(plaintext_bytes, _blk_size)

        # Create AES cipher object and encrypt padded plaintext
        cipher = AES.new(key, AES.MODE_ECB)
        ciphertext_ = cipher.encrypt(padded_plaintext)
        return ciphertext_

    # Store encrypted text in IPFS
    def store_decentralized(self, encrypted_text):
        with ipfshttpclient.connect() as client:
            get_response = client.add_bytes(encrypted_text)
            return get_response

    # Fetch encrypted data from IPFS
    def fetch_data_ipfs(self, ipfs_hash):
        with ipfshttpclient.connect() as client:
            encrypted_text = client.cat(ipfs_hash)
            return encrypted_text

     # Perform AES decryption with custom padding removal
    def perform_aes_decryption(self, key, encrypted_text):
        # Create AES cipher object
        cipher = AES.new(key, AES.MODE_ECB)

        # Decrypt ciphertext
        decrypted_bytes = cipher.decrypt(encrypted_text)

        # Remove custom padding
        decrypted_text = unpad(decrypted_bytes, AES.block_size).decode()
        return decrypted_text

    # Measure the execution time of a function
    def check_effectiveness(self, func, *args, **kwargs):
        _start_time = timeit.default_timer()
        get_result = func(*args, **kwargs)
        _end_time = timeit.default_timer()
        execution_time = _end_time - _start_time
        return get_result, execution_time

    # Run the entire encryption process
    def run_encryption_process(self):
        # Apply DNA encryption to the encryption key
        apply_dna_encryptioned_key, apply_dna_encryptionion_time = self.check_effectiveness(self.apply_dna_encryption)
        print("\nAfter Dna encryption, key:", apply_dna_encryptioned_key)

        # Convert DNA-encrypted key to bytes for AES encryption
        aes_key = apply_dna_encryptioned_key.encode()

        # Perform AES encryption on the plaintext
        encrypted_text, perform_aes_encryptionion_time = self.check_effectiveness(self.perform_aes_encryption, aes_key, self.my_plaintext_)

        # Store encrypted text in IPFS and measure storage time
        ipfs_hash, ipfs_storage_time = self.check_effectiveness(self.store_decentralized, encrypted_text)

        print("\nObtained IPFS Hash:", ipfs_hash)

        # Retrieve encrypted text from IPFS and measure retrieval time
        retrieved_encrypted_text, ipfs_retrieval_time = self.check_effectiveness(self.fetch_data_ipfs, ipfs_hash)
        
        # Decrypt the retrieved text using AES and measure decryption time
        decrypted_text, aes_decryption_time = self.check_effectiveness(self.perform_aes_decryption, aes_key, retrieved_encrypted_text)

        print("\nText after Decryption:", decrypted_text)
        print("\n")
        print("DNA Key Encryption Time:", apply_dna_encryptionion_time, "seconds")
        print("\n")
        print("IPFS Storage Time:", ipfs_storage_time, "seconds")
        print("IPFS Retrieval Time:", ipfs_retrieval_time, "seconds")
        print("\n")
        print("AES Encryption Time:", perform_aes_encryptionion_time, "seconds")
        print("AES Decryption Time:", aes_decryption_time, "seconds")

# Take plaintext input from the user
my_plaintext_ = input("\nEnter the Plaintext :")
sym_key = "pkey"

# Create an instance of the DNAEncryptionSystem class
encryption_system = DNAEncryptionSystem(my_plaintext_, sym_key)

# Run the encryption process
encryption_system.run_encryption_process()
