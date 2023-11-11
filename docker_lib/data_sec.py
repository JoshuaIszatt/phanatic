#!/usr/bin/env python
#

import os
import hashlib

def generate_sha256_hash(path):
    hash_object = hashlib.sha256(path.encode())  # Encodes the string then computes hash
    hex = hash_object.hexdigest()  # Obtains the hexadecimal representation of the hash (hash key)
    return hex

# Finding files