# Import necessary libraries
from Crypto.Cipher import AES,Salsa20
from Crypto.Util.Padding import pad, unpad
from Bio.Seq import Seq
from Crypto.Random import get_random_bytes
import ipfshttpclient
import os
import numpy as np
import binascii
import timeit
import matplotlib.pyplot as plt

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

   # Perform Salsa20 encryption with custom padding
    def perform_salsa20_encryption(self, key, my_plaintext_):
        key = key[:32]  # Salsa20 key size is 256 bits (32 bytes)
        cipher = Salsa20.new(key=key)
        ciphertext_ = cipher.nonce + cipher.encrypt(my_plaintext_.encode())
        return ciphertext_

    # Perform Salsa20 decryption with custom padding removal
    def perform_salsa20_decryption(self, key, encrypted_text):
        key = key[:32]  # Salsa20 key size is 256 bits (32 bytes)
        msg_nonce = encrypted_text[:8]
        ciphertext = encrypted_text[8:]
        cipher = Salsa20.new(key=key, nonce=msg_nonce)
        decrypted_bytes = cipher.decrypt(ciphertext)
        return decrypted_bytes

    # Measure the execution time of a function
    def check_effectiveness(self, func, *args, **kwargs):
        _start_time = timeit.default_timer()
        get_result = func(*args, **kwargs)
        _end_time = timeit.default_timer()
        execution_time = _end_time - _start_time
        return get_result, execution_time

    # Run the entire encryption process with AES
    def run_encryption_process(self):
        # Apply DNA encryption to the encryption key
        apply_dna_encryptioned_key, apply_dna_encryptionion_time = self.check_effectiveness(self.apply_dna_encryption)
        print("\nAfter DNA encryption, key:", apply_dna_encryptioned_key)

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

        return (
            perform_aes_encryptionion_time,
            aes_decryption_time
        )

    # Run the entire encryption process with salsa20
    def run_encryption_process_with_salsa20(self):
        # Apply DNA encryption to the encryption key
        apply_dna_encryptioned_key, apply_dna_encryptionion_time = self.check_effectiveness(self.apply_dna_encryption)
        print("\nAfter DNA encryption, key:", apply_dna_encryptioned_key)

        # Convert DNA-encrypted key to bytes for salsa20 encryption
        salsa20_key = apply_dna_encryptioned_key.encode()

        # Perform salsa20 encryption on the plaintext
        encrypted_text, perform_salsa20_encryption_time = self.check_effectiveness(self.perform_salsa20_encryption, salsa20_key, self.my_plaintext_)

        # Store encrypted text in IPFS and measure storage time
        ipfs_hash, ipfs_storage_time = self.check_effectiveness(self.store_decentralized, encrypted_text)

        print("\nObtained IPFS Hash:", ipfs_hash)

        # Retrieve encrypted text from IPFS and measure retrieval time
        retrieved_encrypted_text, ipfs_retrieval_time = self.check_effectiveness(self.fetch_data_ipfs, ipfs_hash)

        # Decrypt the retrieved text using salsa20 and measure decryption time
        decrypted_text, salsa20_decryption_time = self.check_effectiveness(self.perform_salsa20_decryption, salsa20_key, retrieved_encrypted_text)

        print("\nText after Decryption:", decrypted_text)
        print("\n")
        print("DNA Key Encryption Time:", apply_dna_encryptionion_time, "seconds")
        print("\n")
        print("IPFS Storage Time:", ipfs_storage_time, "seconds")
        print("IPFS Retrieval Time:", ipfs_retrieval_time, "seconds")
        print("\n")
        print("Salsa20 Encryption Time:", perform_salsa20_encryption_time, "seconds")
        print("Salsa20 Decryption Time:", salsa20_decryption_time, "seconds")

        return (
            perform_salsa20_encryption_time,
            salsa20_decryption_time
        )

# Take plaintext input from the user
my_plaintext_ = input("\nEnter the Plaintext :")
sym_key = "pkey"

# Create an instance of the DNAEncryptionSystem class
encryption_system = DNAEncryptionSystem(my_plaintext_, sym_key)

# Run the encryption process and get the times for AES
aes_times = encryption_system.run_encryption_process()

# Run the encryption process with salsa20 and get the times
salsa20_times = encryption_system.run_encryption_process_with_salsa20()

# Plotting the comparison graph
labels = ['Encryption', 'Decryption']

# Extracting times for AES and Salsa20
aes_times = list(aes_times)
salsa20_times = list(salsa20_times)

x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, aes_times, width, label='AES')
rects2 = ax.bar(x + width/2, salsa20_times, width, label='Salsa20')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Time (seconds)')
ax.set_title('Comparison of AES and Salsa20')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

# Add values on top of the bars
def add_values(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(round(height, 3)),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

add_values(rects1)
add_values(rects2)

fig.tight_layout()

plt.show()