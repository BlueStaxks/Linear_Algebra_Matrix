The code is designed to detect fake shares and successfully reconstruct the secret key K.

This detection algorithm is very robust. It can reconstruct the key even with 11 valid and 1000 fake shares. (t=10, s=1001, prime=524287)

This method requires a (s+1) by (m) matrix where s is maximum number of fake shares and m is number of points including valid and fake shares.

Reciever uses the matrix to find valid shares and only use them to reconstruct K.

If number of real shares are less than t+1, the system cannot reconstruct K because the matrix does not have information of y coordinates so it is secure. 


*Beware, that this method is not guaranteed to succeed, but highly probable to succeed.

(mathematical explanation will be updated)