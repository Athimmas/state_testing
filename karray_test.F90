program karray
implicit none

integer ip_array(2,2,2),op_array(4,2)
integer i,j,k

do i=1,2
ip_array(:,:,i)=i
enddo

op_array=reshape(ip_array,(/4,2/))

print *,ip_array
print *,op_array

end program karray
