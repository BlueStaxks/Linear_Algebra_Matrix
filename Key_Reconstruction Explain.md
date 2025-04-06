The code is designed to detect fake shares and successfully reconstruct the secret key K.

This detection algorithm is very robust. It can reconstruct the key even with 31 valid and 900 fake shares.

This method requires a (s+1) by (m) matrix where s is maximum number of fake shares and m is number of points including valid and fake shares.

Reciever uses the matrix to find valid shares and only use them to reconstruct K.
