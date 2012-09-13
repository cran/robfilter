#ifndef CIRCULARRRAY_H
#define CIRCULARRRAY_H

#include <stdexcept>
#include <iostream>

template<typename T>
class CircularArray {
      int capacity;
      
      T* data;
      
      /**
       * Where the next element should be inserted.
       */
      int next;
      
      /**
       * The number of elements that have been inserted.
       */
      int len;
      
      //Hide the copy constructor and assignment operator
      CircularArray(const CircularArray&);
      CircularArray& operator=(CircularArray&);
      
      /**
       * Index of the oldest element in the array.
       */
      int indexOfOldest() const {
          //If len > 0, next points one index past the newest element % capacity.
          //Thus we must go back len places % capacity to reach the oldest, eg,
          //if len == 1, we must count back 1 place.
          //-1 == capacity - 1 (mod capacity). 
          //Thus, if -capacity <= next - len < 0, we will find the oldest element
          //at capacity + next - len
          int oldest = next - len;
          if (oldest < 0) {
             oldest += capacity;
          }
          return oldest;
      }
      
      public:
      CircularArray(int capacity);
      ~CircularArray();
      
      /**
       * Add an element to the "front" of the array.
       * If size() == capacity already,
       * remove the oldest element from the array.
       */
      void append(const T& x);
      
      /**
       * Return the ith element in this array, the elements
       * being arranged in order starting with the oldest at 0.
       * Throws an out_of_range exception if there is no such element.
       */
      T& get(int i);
      
      /**
       * Return the value of the oldest element
       * in the array.
       * Throw an out_of_range exception if there
       * is no such element.
       */
      T oldest() const;
      
      /**
       * Remove the least recently inserted element.
       * Do nothing if there is no such element.
       */
      void removeOldest();
      
      /**
       * Gets the number of elements currently in the array.
       * This is always <= capacity.
       */
      int size() const;
      
      /**
       * Copies the elements currently in the array into sink
       * in the order in which they were inserted.
       * Does nothing if size() == 0.
       */
      void copyInto(T* sink) const;
};

template<typename T>
CircularArray<T>::CircularArray(int maxSize) : capacity(maxSize), next(0), len(0)
{
 if (capacity > 0) {
   data = new T[capacity];
 }
 else if (capacity == 0) {
   //Ensuring that data always has room for at least
   //one element allows us to safely handle "appending"
   //an element to an array with capacity zero
   //without burdening the append method with extra
   //code.
   data = new T[1];
 }
 else {
   throw std::invalid_argument("capacity < 0");
 }
}

template<typename T>
CircularArray<T>::~CircularArray<T>()
{
 delete[] data;
}

/**
 * Add an element to the "front" of the array.
 * If size() == capacity already,
 * remove the oldest element from the array.
 */
template<typename T>
void CircularArray<T>::append(const T& x)
{
 data[next] = x;
 ++next;
 if (next >= capacity) {
   next = 0;
 }
 if (len < capacity) {
    ++len;
 }
}

/**
 * Return the ith element in this array, the elements
 * being arranged in order starting with the oldest at 0.
 * Throws an out_of_range exception if there is no such element.
 */
template<typename T>
T& CircularArray<T>::get(int i)
{
  //XXX! This seems pretty heavy compared to a simple array access.
  //Should we throw out the range check?
  if (i < 0 || i >= len) {
    throw std::out_of_range("index out of range");
  }
  else {
    const int oldest = indexOfOldest();
    return data[(oldest + i) % capacity];
  }
}

/**
 * Remove the least recently inserted element.
 * Do nothing if there is no such point.
 */
template<typename T>
void CircularArray<T>::removeOldest()
{
 if (len > 0) {
    --len;
 }
}

/**
 * Return the value of the oldest element
 * in the array.
 * Throw an out_of_range exception if there
 * is no such element.
 */
template<typename T>
T CircularArray<T>::oldest() const
{
  //std::cout << "entering oldest(),   len=" << len << "\n";
            	
  if (len > 0) {
   //  std::cout << "leaving oldest() --> normal way \n";
     return data[indexOfOldest()];
  }
  else {
 //   std::cout << "leaving oldest() --> throwing exception \n";
    throw std::out_of_range("no such element");
  }
}

/**
 * Gets the number of elements currently in the array.
 * This is always <= capacity.
 */
template<typename T>
int CircularArray<T>::size() const
{
    return len;
}

/**
 * Copies the elements currently in the array into sink
 * in the order in which they were inserted.
 * Does nothing if size() == 0.
 */
template<typename T>
void CircularArray<T>::copyInto(T* sink) const
{
     //Note: Could use memcpy here.
     
     //Index of the oldest element not yet copied.
     int oldest = indexOfOldest();
     for (int i = 0; i != len; ++i, ++oldest) {
         sink[i] = data[oldest % capacity];
     }
}

#endif
