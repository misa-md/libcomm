//
// Created by genshen on 2019-03-16.
//

#ifndef COMM_PACKER_H
#define COMM_PACKER_H

namespace comm {
  /**
   * \tparam RT type of region type
   * \brief abstract class for data packer.
   */
  template <typename T> class Packer {
  public:
    typedef T pack_date_type;

    /**
     * count the length to be send to neighbour process.
     * \param dimension 0,1,2. the id of dimensions.
     * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
     * \return
     */
    virtual const unsigned long sendLength(const int dimension, const int direction) = 0;

    /**
     * This function will be called before sending data.
     * \param buffer the empty buffer to be send to its neighbor process, you can fill the buffer here.
     * \param send_size the length to be send.
     * \param dimension 0,1,2. the id of dimensions.
     * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
     */
    virtual void onSend(T buffer[], const unsigned long send_len, const int dimension, const int direction) = 0;

    /**
     * This function will be called after the receive finished.
     * We can unpack data here.
     * \param buffer the buffer of received data.
     * \param receive_len the length of received data.
     * \param dimension 0,1,2. the id of dimensions.
     * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
     */
    virtual void onReceive(T buffer[], const unsigned long receive_len, const int dimension, const int direction) = 0;

    /**
     * \brief this function will be called after all communication finished.
     */
    virtual void onFinish(){};
  };
} // namespace comm

#endif // COMM_PACKER_H
