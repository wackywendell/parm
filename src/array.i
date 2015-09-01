%{
#include <boost/array.hpp>
%}

namespace boost {

template<typename T, std::size_t N>
class array {
public:
    array();
    ~array();

    void swap(array& b);

    %extend {
        array(){
            return new array<T, N>();
        }
        
        T __getitem__(unsigned int idx){
            if(idx >= N){
                PyErr_SetString(PyExc_IndexError, "Index invalid.");
                return (*self)[0];
            }
            return (*self)[idx];
        }
        
        void __setitem__(unsigned int idx, T val){
            (*self)[idx]=val;
        }
        
        unsigned int __len__(){ return N;};
        
        ~array() {
            %delete(self);
        }
        
        %insert("python") %{
            def __iter__(self):
                for i in range(len(self)):
                    yield self[i]
            
            def __str__(self):
                return str(tuple(self))
        
            def __repr__(self):
                return repr(tuple(self))
        %}
    };
};

}
