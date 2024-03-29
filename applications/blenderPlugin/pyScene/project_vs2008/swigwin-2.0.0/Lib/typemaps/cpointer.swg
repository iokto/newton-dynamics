/* -----------------------------------------------------------------------------
 * cpointer.swg
 *
 * This library file contains macros that can be used to manipulate simple
 * pointer objects.
 *
 * ----------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
 * %pointer_class(type,name)
 *
 * Places a simple proxy around a simple type like 'int', 'float', or whatever.
 * The proxy provides this interface:
 *
 *       class type {
 *       public:
 *           type();
 *          ~type();
 *           type value();
 *           void assign(type value);
 *       };
 *         
 * Example:
 *
 *    %pointer_class(int, intp);
 *
 *    int add(int *x, int *y) { return *x + *y; }
 *
 * In python (with proxies)
 *
 *    >>> a = intp()
 *    >>> a.assign(10)
 *    >>> a.value()
 *    10
 *    >>> b = intp()
 *    >>> b.assign(20)
 *    >>> print add(a,b)
 *    30
 *
 * As a general rule, this macro should not be used on class/structures that
 * are already defined in the interface.
 * ----------------------------------------------------------------------------- */


%define %pointer_class(TYPE, NAME)
%{
typedef TYPE NAME;
%}

typedef struct {
} NAME;

%extend NAME {
  NAME() {
    return %new_instance(TYPE);
  }
  ~NAME() {
    if (self) %delete(self);
  }
}

%extend NAME {

  void assign(TYPE value) {
    *self = value;
  }
  TYPE value() {
    return *self;
  }
  TYPE * cast() {
    return self;
  }
  static NAME * frompointer(TYPE *t) {
    return (NAME *) t;
  }
}

%types(NAME = TYPE);

%enddef

/* ----------------------------------------------------------------------------- 
 * %pointer_functions(type,name)
 *
 * Create functions for allocating/deallocating pointers.   This can be used
 * if you don't want to create a proxy class or if the pointer is complex.
 *
 *    %pointer_functions(int, intp)
 *
 *    int add(int *x, int *y) { return *x + *y; }
 *
 * In python (with proxies)
 *
 *    >>> a = copy_intp(10)
 *    >>> intp_value(a)
 *    10
 *    >>> b = new_intp()
 *    >>> intp_assign(b,20)
 *    >>> print add(a,b)
 *    30
 *    >>> delete_intp(a)
 *    >>> delete_intp(b)
 * 
 * ----------------------------------------------------------------------------- */

%define %pointer_functions(TYPE,NAME)
%{
  static TYPE *new_##NAME() { 
    return %new_instance(TYPE);
  }
  
  static TYPE *copy_##NAME(TYPE value) { 
    return %new_copy(value, TYPE);
  }

  static void delete_##NAME(TYPE *self) { 
    if (self) %delete(self);
  }

  static void NAME ##_assign(TYPE *self, TYPE value) {
    *self = value;
  }

  static TYPE NAME ##_value(TYPE *self) {
    return *self;
  }
%}

TYPE *new_##NAME();
TYPE *copy_##NAME(TYPE value);
void  delete_##NAME(TYPE *self);
void  NAME##_assign(TYPE *self, TYPE value);
TYPE  NAME##_value(TYPE *self);

%enddef

/* -----------------------------------------------------------------------------
 * %pointer_cast(type1,type2,name)
 *
 * Generates a pointer casting function.
 * ----------------------------------------------------------------------------- */

%define %pointer_cast(TYPE1,TYPE2,NAME)
%inline %{
TYPE2 NAME(TYPE1 x) {
   return %static_cast(x, TYPE2);
}
%}
%enddef








