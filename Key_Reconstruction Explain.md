The code is designed to detect fake shares and successfully reconstruct the secret key K.

This detection algorithm is very robust. It can reconstruct the key even with 11 valid and 1000 fake shares. (t=10, s=1001, prime=524287)

This method requires a (s+1) by (m) matrix where s is maximum number of fake shares and m is number of points including valid and fake shares.

Reciever uses the matrix to find valid shares and only use them to reconstruct K.
