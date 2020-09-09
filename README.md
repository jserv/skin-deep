# skin-deep: Generate faces with smoother skin

This program takes effective steps to generate pixels with smooth skin:
1. Skin detection: detect image pixels and regions that contain skin-tone color.
2. Edge detection: obtain an edge map for applying noise reduction over skin.

## Build and Verify

Run `make` and you should see the generated `main` file. Its usage:
```
./main <infile.jpg>
```

The output filename is `out.jpg`.

## License

`skin-deep`is released under the BSD 2 clause license. Use of this source code is governed by
a BSD-style license that can be found in the LICENSE file.

External source code:
* `stb_image.h` and `stb_image_write.h`: Taken from [stb](https://github.com/nothings/stb),
  MIT License or public domain.
