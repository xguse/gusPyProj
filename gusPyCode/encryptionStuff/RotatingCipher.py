class RotatingCipher(object):
  def __init__(self, alphabet, cipher):
    self._alphabet = alphabet
    self._cipher = cipher

  # to be overridden by subclasses for different rotation
  # schemes.
  def _RotateCipher(self, c):
    pass

  def CryptNextChar(self, c):
    """Encrypt a single character and rotate the cipher."""
    # find position of "c" in the alphabet.
    index = self._alphabet.find(c)
    # encrypted character is that position plus the rotation,
    # modulo the size of the alphabet.    
    crypt = self._cipher[(index + self._rot) % len(self._alphabet)]
    self._RotateCipher(c)
    return crypt

  def Encrypt(self, s):
    """Encrypt a message string."""
    self._rot = 0
    result = ""
    for c in s:
      result += self.CryptNextChar(c)
    return result

  def DecryptNextChar(self, c):
    """Decrypt a single character and rotate the cipher."""
    # to decrypt of character, find index of the crypted character in
    # the cipher, and subtract the rotation from it, giving the
    # index of the decrypted character in the alphabet.
    index = self._cipher.find(c)
    index -= self._rot
    if index < 0: index = len(self._alphabet) + index
    self._RotateCipher(self._alphabet[index])
    return self._alphabet[index]

  def Decrypt(self, s):
    """Decrypt a string encrypted by this cipher."""
    self._rot = 0
    result = ""
    for c in s:
      result += self.DecryptNextChar(c)
    return result


class FixedRotatingCipher(RotatingCipher):
  """An implementation of a cipher using a rotation strategy 
     that shifts the cipher by a fixed increment after each
     character."""
  def __init__(self, alphabet, cipher, increment):
    super(FixedRotatingCipher, self).__init__(alphabet, cipher)
    self._rot = 0
    self._increment = increment

  def _RotateCipher(self, c):
    self._rot = (self._rot + self._increment) % len(self._alphabet)



class VariableRotatingCipher(RotatingCipher):
    """An implementation of a cipher using a function of the
       last character as the rotator. Each rotation shifts by
       the character index times three."""
    def __init__(self, alphabet, cipher):
        super(VariableRotatingCipher, self).__init__(alphabet, cipher)

    def _RotateCipher(self, c):
        index = self._alphabet.find(c)
        self._rot = (self._rot + (3*index)) % len(self._alphabet)


fcipher = FixedRotatingCipher(
     "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.,- ", 
     "KLMNijklSTUVabcduvwxyzABCDEFefghGHIJ.,- mnopOPQRqrstWXYZ", 3)

en = fcipher.Encrypt("Mark has a big nose.")
de = fcipher.Decrypt(en)

print "Encrypted = %s" % en
print "Decrypted = %s" % de

vcipher = VariableRotatingCipher(
     "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ.,- ", 
     "KLMNijklSTUVabcduvwxyzABCDEFefghGHIJ.,- mnopOPQRqrstWXYZ")

ven = vcipher.Encrypt("Mark has a big nose.")
vde = vcipher.Decrypt(ven)

print "Variable encrypted = %s" % ven
print "Variable decrypted = %s" % vde