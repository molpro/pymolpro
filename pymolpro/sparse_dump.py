import struct
import numpy as np


def sparse_dump_get(filename, offset=0):
    """
    Unpack a sparse_dump representation on file

    :param filename: The filename of the sparse_dump
    :type filename: str
    :param offset: Offset of the sparse_dump from the beginning of the file
    :type offset: int
    :return: uncompressed values, length of sparse_dump
    :rtype: np.array, int
    """
    size_t_length = 8
    type_length_t_length = 1
    size_t_format = '<q'
    type_length_t_format = '<b'

    with open(filename, 'rb') as f:
        f.seek(offset)
        assert struct.calcsize(size_t_format) == size_t_length
        assert struct.calcsize(type_length_t_format) == type_length_t_length
        full_length, = struct.unpack(size_t_format, f.read(size_t_length))
        segment_skip_size, = struct.unpack(type_length_t_format, f.read(type_length_t_length))
        for sparse_int_type in [np.uint8, np.uint16, np.uint32, np.uint64]:
            if segment_skip_size == sparse_int_type(0).itemsize: break
        value_size, = struct.unpack(type_length_t_format, f.read(type_length_t_length))
        for sparse_float_type in [np.float32, np.float64]:
            if value_size == sparse_float_type(0).itemsize: break
        number_of_segments, = struct.unpack(size_t_format, f.read(size_t_length))
        number_of_values, = struct.unpack(size_t_format, f.read(size_t_length))
        segment_skip = np.fromfile(f, count=number_of_segments, dtype=sparse_int_type)
        segment_length = np.fromfile(f, count=number_of_segments, dtype=sparse_int_type)
        values = np.fromfile(f, count=number_of_values, dtype=sparse_float_type)
        length = f.tell()-offset
        buffer = np.zeros(full_length, dtype=np.float64)
        offset = 0
        packed_offset = 0
        for segment in range(number_of_segments):
            offset += segment_skip[segment]
            buffer[offset:offset + segment_length[segment]] = values[
                                                              packed_offset:packed_offset + segment_length[segment]]
            packed_offset += segment_length[segment]
    return buffer, length


if __name__ == "__main__":
    import pathlib
    file = pathlib.Path(__file__).parent.parent / 'example.xml-sidecar'
    buffer, length = sparse_dump_get(file)
    print(length)
    print(buffer)
    offset = length
    buffer, length = sparse_dump_get(file, offset=offset)
    print(length)
    print(buffer)
    print(length)
    print(buffer)
