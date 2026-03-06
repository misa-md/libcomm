//
// Created by genshen on 2019-03-16.
//

#ifndef COMM_PACKER_H
#define COMM_PACKER_H

#include <vector>

namespace comm {
  /**
   * perform the forward communication with std::vector as packer container (without the sendLength api).
   * @tparam T the datatype used for sending and receiving array.
   */
  template <typename T> class Packer {
  public:
    using pack_date_type = T;

    virtual ~Packer() = default;

    /**
     * this function is called before performing communication.
     */
    virtual void initialize() = 0;

    /**
     * pack the data into the @param data for MPI sending on dimension @param dimension and direction @param direction
     * @param data the empty data buffer for the sending data.
     * @param dimension the dimension of x(0) ,y (1) , z(2)
     * @param direction the sending direction: lower or higher.
     * @return
     */
    virtual const std::size_t pack(std::vector<T> data, const int dimension, const int direction) = 0;

    /**
     * unpack data after performing MPI_Recv (e.g. save to another temp buffer).
     * you can implement the unpacking details in the function.
     *
     * @param data the received data from MPI.
     * @param data_len the length of the received data.
     * @param dimension the dimension of x(0) ,y (1) , z(2).
     * @param direction the sending direction: lower or higher.
     */
    virtual const void unpack(const std::vector<T> data, const std::size_t data_len, const int dimension,
                              const int direction) = 0;

    /**
     * this function is called when the communication for all dimension and all direction is finished.
     */
    virtual void done() = 0;
  };

  /**
   * \tparam T type of region type
   * \brief abstract class for data packer.
   *
   * It is similar as \class Packer, but it uses the native C-array style interface and calculate the data length via
   * \memberof sendLength before packing data.
   */
  template <typename T> class NativePacker : Packer<T> {
  protected:
    void initialize() override {}

    inline const std::size_t pack(std::vector<T> data, const int dimension, const int direction) override {
      std::size_t send_len = sendLength(dimension, direction);
      data.resize(send_len); // malloc storage
      onSend(data.data(), send_len, dimension, direction);
      return send_len;
    }

    inline const void unpack(const std::vector<T> data, const std::size_t data_len, const int dimension,
                             const int direction) override {
      onReceive(data.data(), data_len, dimension, direction);
    }

    inline void done() override { onFinish(); }

  public:
    using pack_date_type = T;

    /**
     * count the length to be sent to neighbour process.
     * \param dimension 0,1,2. the id of dimensions.
     * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
     * \return
     */
    virtual const std::size_t sendLength(const int dimension, const int direction) = 0;

    /**
     * This function will be called before sending data.
     * \param buffer the empty buffer to be sent to its neighbor process, you can fill the buffer here.
     * \param send_len the length to be sent.
     * \param dimension 0,1,2. the id of dimensions.
     * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
     */
    virtual void onSend(T buffer[], const std::size_t send_len, const int dimension, const int direction) = 0;

    /**
     * This function will be called after the received finished.
     * We can unpack data here.
     * \param buffer the buffer of received data.
     * \param receive_len the length of received data.
     * \param dimension 0,1,2. the id of dimensions.
     * \param direction DIR_LOWER or DIR_HIGHER, the direction id.
     */
    virtual void onReceive(T buffer[], const std::size_t receive_len, const int dimension, const int direction) = 0;

    /**
     * \brief this function will be called after all communication finished.
     */
    virtual void onFinish() {};
  };
} // namespace comm

#endif // COMM_PACKER_H
