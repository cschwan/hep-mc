#include "gtest/gtest.h"

#ifndef HEP_USE_MPI
#include "hep/mc.hpp"
#else
#include "hep/mc-mpi.hpp"
#endif

#include <cmath>
#include <cstddef>
//#include <cstdio>
#include <tuple>
#include <vector>

template <typename T>
T arctan(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
	T const x = T(10.0) * point.point().at(0) - T(5.0);
	T const f = std::atan(x);

	projector.add(0, x, f);

	return f;
}

template <typename T>
std::vector<std::tuple<T, T>> reference_results();

template <>
std::vector<std::tuple<float, float>> reference_results() {
	return {
		std::make_tuple(-0x9.1a251p-6f, 0x8.ff77ep-11f),
		std::make_tuple(-0x8.dd8a8p-6f, 0x8.de413p-11f),
		std::make_tuple(-0x8.8cdcfp-6f, 0x8.b27b8p-11f),
		std::make_tuple(-0x8.52c7ap-6f, 0x8.918bfp-11f),
		std::make_tuple(-0x8.44fc8p-6f, 0x8.86ce0p-11f),
		std::make_tuple(-0x9.4b79ap-6f, 0x9.04c1ep-11f),
		std::make_tuple(-0x8.89d59p-6f, 0x8.a1c78p-11f),
		std::make_tuple(-0x8.d4c2fp-6f, 0x8.c28a6p-11f),
		std::make_tuple(-0x8.b7e7fp-6f, 0x8.afc15p-11f),
		std::make_tuple(-0x9.08179p-6f, 0x8.d2467p-11f),
		std::make_tuple(-0x8.774eap-6f, 0x8.862e2p-11f),
		std::make_tuple(-0x8.66ca5p-6f, 0x8.78d6bp-11f),
		std::make_tuple(-0x8.75c59p-6f, 0x8.7aca0p-11f),
		std::make_tuple(-0x8.1bd69p-6f, 0x8.48460p-11f),
		std::make_tuple(-0x8.6dba1p-6f, 0x8.6b2f9p-11f),
		std::make_tuple(-0x8.151f2p-6f, 0x8.38af7p-11f),
		std::make_tuple(-0x8.7146fp-6f, 0x8.5fdc7p-11f),
		std::make_tuple(-0xf.80e5dp-7f, 0x8.00515p-11f),
		std::make_tuple(-0x8.01430p-6f, 0x8.19d39p-11f),
		std::make_tuple(-0x8.0edfcp-6f, 0x8.18d95p-11f),
		std::make_tuple(-0xf.f000dp-7f, 0x8.04f7fp-11f),
		std::make_tuple(-0x8.1b0d9p-6f, 0x8.0d5e9p-11f),
		std::make_tuple(-0xf.6014ap-7f, 0xf.9d918p-12f),
		std::make_tuple(-0xf.787f0p-7f, 0xf.95b46p-12f),
		std::make_tuple(-0xf.90c24p-7f, 0xf.8c153p-12f),
		std::make_tuple(-0xf.e621bp-7f, 0xf.9e4f9p-12f),
		std::make_tuple(-0xe.ca66dp-7f, 0xe.f9c49p-12f),
		std::make_tuple(-0xf.011d4p-7f, 0xe.fab66p-12f),
		std::make_tuple(-0xe.15acap-7f, 0xe.68a50p-12f),
		std::make_tuple(-0xe.4ae5ep-7f, 0xe.640cbp-12f),
		std::make_tuple(-0xd.ffad5p-7f, 0xe.1dcdcp-12f),
		std::make_tuple(-0xd.8e405p-7f, 0xd.c0e54p-12f),
		std::make_tuple(-0xd.eb75cp-7f, 0xd.c7e70p-12f),
		std::make_tuple(-0xd.2942fp-7f, 0xd.3d0e0p-12f),
		std::make_tuple(-0xc.55e94p-7f, 0xc.a3ba1p-12f),
		std::make_tuple(-0xc.76732p-7f, 0xc.8132dp-12f),
		std::make_tuple(-0xb.f7412p-7f, 0xc.09595p-12f),
		std::make_tuple(-0xc.24a73p-7f, 0xb.e0e24p-12f),
		std::make_tuple(-0xb.661d3p-7f, 0xb.3f6a2p-12f),
		std::make_tuple(-0xa.7412cp-7f, 0xa.7aaa8p-12f),
		std::make_tuple(-0x9.9f8c3p-7f, 0x9.bed6ap-12f),
		std::make_tuple(-0x8.86c26p-7f, 0x8.d4d6dp-12f),
		std::make_tuple(-0xf.cffe7p-8f, 0x8.224f5p-12f),
		std::make_tuple(-0xe.f96a8p-8f, 0xe.f937ep-13f),
		std::make_tuple(-0xd.0aaa0p-8f, 0xd.122e2p-13f),
		std::make_tuple(-0xa.dbeeep-8f, 0xa.ef6d5p-13f),
		std::make_tuple(-0x8.87151p-8f, 0x8.a70e3p-13f),
		std::make_tuple(-0xb.c7e61p-9f, 0xc.5366bp-14f),
		std::make_tuple(-0xe.7f454p-10f, 0xf.2b522p-15f),
		std::make_tuple(-0x9.e81a1p-11f, 0xb.b6ac3p-16f),
		std::make_tuple( 0xb.0d473p-11f, 0xc.6da31p-16f),
		std::make_tuple( 0xf.1d2bdp-10f, 0xf.90a26p-15f),
		std::make_tuple( 0xc.36107p-9f, 0xc.84340p-14f),
		std::make_tuple( 0x8.cd1dap-8f, 0x8.c5074p-13f),
		std::make_tuple( 0xb.3700ep-8f, 0xb.1ab53p-13f),
		std::make_tuple( 0xc.c0cb8p-8f, 0xc.e9a0fp-13f),
		std::make_tuple( 0xe.51189p-8f, 0xe.a3c14p-13f),
		std::make_tuple( 0x8.2719ap-7f, 0x8.3fe5cp-12f),
		std::make_tuple( 0x8.e9f49p-7f, 0x9.08517p-12f),
		std::make_tuple( 0x9.e8637p-7f, 0x9.e361dp-12f),
		std::make_tuple( 0xa.ac7c0p-7f, 0xa.96da7p-12f),
		std::make_tuple( 0xb.5c74dp-7f, 0xb.3a0d8p-12f),
		std::make_tuple( 0xb.26730p-7f, 0xb.63c60p-12f),
		std::make_tuple( 0xc.64654p-7f, 0xc.3ebf5p-12f),
		std::make_tuple( 0xc.4a578p-7f, 0xc.6b570p-12f),
		std::make_tuple( 0xd.256a1p-7f, 0xd.0a68bp-12f),
		std::make_tuple( 0xd.3cf18p-7f, 0xd.4696bp-12f),
		std::make_tuple( 0xd.b0546p-7f, 0xd.aa8a3p-12f),
		std::make_tuple( 0xe.64bb7p-7f, 0xe.2abf3p-12f),
		std::make_tuple( 0xd.ac3b8p-7f, 0xd.f37f6p-12f),
		std::make_tuple( 0xe.27516p-7f, 0xe.52ca1p-12f),
		std::make_tuple( 0xe.09d88p-7f, 0xe.624e6p-12f),
		std::make_tuple( 0xe.40524p-7f, 0xe.9a1dcp-12f),
		std::make_tuple( 0xe.80f8bp-7f, 0xe.d46bcp-12f),
		std::make_tuple( 0xe.d48c9p-7f, 0xf.17cfep-12f),
		std::make_tuple( 0xf.2b099p-7f, 0xf.5998ep-12f),
		std::make_tuple( 0xe.e6bcdp-7f, 0xf.4cce9p-12f),
		std::make_tuple( 0x8.287a7p-6f, 0x8.0a911p-11f),
		std::make_tuple( 0xf.e53aap-7f, 0xf.f2c11p-12f),
		std::make_tuple( 0x8.108c8p-6f, 0x8.1133fp-11f),
		std::make_tuple( 0x8.1ce4ap-6f, 0x8.1fa05p-11f),
		std::make_tuple( 0xf.f164dp-7f, 0x8.1562dp-11f),
		std::make_tuple( 0xf.95c79p-7f, 0x8.05a81p-11f),
		std::make_tuple( 0x8.59d35p-6f, 0x8.5421ap-11f),
		std::make_tuple( 0x8.3904fp-6f, 0x8.4aaf6p-11f),
		std::make_tuple( 0x8.436e8p-6f, 0x8.56440p-11f),
		std::make_tuple( 0x8.7101ap-6f, 0x8.72c44p-11f),
		std::make_tuple( 0x8.55d9bp-6f, 0x8.6b08cp-11f),
		std::make_tuple( 0x8.48a01p-6f, 0x8.69bf4p-11f),
		std::make_tuple( 0x8.82347p-6f, 0x8.8ba27p-11f),
		std::make_tuple( 0x8.e30a5p-6f, 0x8.c047bp-11f),
		std::make_tuple( 0x8.a2300p-6f, 0x8.a517ap-11f),
		std::make_tuple( 0x8.ad4fdp-6f, 0x8.af1efp-11f),
		std::make_tuple( 0x8.7ebd0p-6f, 0x8.9c2dep-11f),
		std::make_tuple( 0x8.1edbcp-6f, 0x8.6f8c1p-11f),
		std::make_tuple( 0x8.8c0cap-6f, 0x8.aace2p-11f),
		std::make_tuple( 0x8.5ddb6p-6f, 0x8.97285p-11f),
		std::make_tuple( 0x8.bded0p-6f, 0x8.cb028p-11f),
		std::make_tuple( 0x8.a5900p-6f, 0x8.c261cp-11f),
		std::make_tuple( 0x8.eade3p-6f, 0x8.e82d9p-11f)
	};
}

template <>
std::vector<std::tuple<double, double>> reference_results() {
	return {
		std::make_tuple(-0x9.591830051bfd8p-6, 0x9.1e0fc505ad23p-11),
		std::make_tuple(-0x8.db45b19d91c18p-6, 0x8.dd1e17cb0aaap-11),
		std::make_tuple(-0x8.c4cb5a7710dp-6, 0x8.ce82fd5e6869p-11),
		std::make_tuple(-0x8.b93b2e9f26ac8p-6, 0x8.c528020814848p-11),
		std::make_tuple(-0x8.832c4ddd456e8p-6, 0x8.a6587d58a8468p-11),
		std::make_tuple(-0x9.06de7a096282p-6, 0x8.e38f934f3f5f8p-11),
		std::make_tuple(-0x8.e65ae9841ae68p-6, 0x8.cf923e7a71128p-11),
		std::make_tuple(-0x8.9ded366d7106p-6, 0x8.a779e8890e898p-11),
		std::make_tuple(-0x8.af403caa0e6a8p-6, 0x8.ab875af7e015p-11),
		std::make_tuple(-0x8.dc8dff53f7c1p-6, 0x8.bd2416288b5cp-11),
		std::make_tuple(-0x8.f2f5199298ad8p-6, 0x8.c2f9cd31a5568p-11),
		std::make_tuple(-0x8.8b67a69f4c8c8p-6, 0x8.8b012260a5a48p-11),
		std::make_tuple(-0x8.2ab1d3bf0ed28p-6, 0x8.553bbcae42aap-11),
		std::make_tuple(-0x8.1985f7077c45p-6, 0x8.4702beccff96p-11),
		std::make_tuple(-0x8.60e183556c488p-6, 0x8.64c65610146b8p-11),
		std::make_tuple(-0xf.a320969236cb8p-7, 0x8.1660ca7f2e09p-11),
		std::make_tuple(-0x8.6cebbe36bcedp-6, 0x8.5da41cb7ee1dp-11),
		std::make_tuple(-0xf.fa1ce4eec1488p-7, 0x8.1f2303db1b2f8p-11),
		std::make_tuple(-0x8.0770561d6f87p-6, 0x8.1ce6d8a0a635p-11),
		std::make_tuple(-0xf.a22af789693ep-7, 0xf.f39ac767f76bp-12),
		std::make_tuple(-0xf.c32535fa1e688p-7, 0xf.f384cf8af77c8p-12),
		std::make_tuple(-0xf.edc2ef78ce73p-7, 0xf.f73786634669p-12),
		std::make_tuple(-0xf.b8eba35bbef4p-7, 0xf.ca5b8d8693cap-12),
		std::make_tuple(-0xf.884ebeb2366b8p-7, 0xf.9d8e6812e7368p-12),
		std::make_tuple(-0xf.fb33b8f5de32p-7, 0xf.c09cd0c9d45bp-12),
		std::make_tuple(-0xf.e67081a283048p-7, 0xf.9e9ca3b517dap-12),
		std::make_tuple(-0xe.38cdf0bd20a98p-7, 0xe.affedea00ec48p-12),
		std::make_tuple(-0xf.3d207535a71e8p-7, 0xf.1808c08f37c28p-12),
		std::make_tuple(-0xe.c075c18160c8p-7, 0xe.bded27e11fb88p-12),
		std::make_tuple(-0xe.140b435cedc28p-7, 0xe.489f241e1b6d8p-12),
		std::make_tuple(-0xe.9d68ad090418p-7, 0xe.6b717b25e7d28p-12),
		std::make_tuple(-0xd.87678ee3e24dp-7, 0xd.bd8edc23bde3p-12),
		std::make_tuple(-0xd.837eeaccafc18p-7, 0xd.944cd80a5a64p-12),
		std::make_tuple(-0xc.9110a1478ef68p-7, 0xc.efec35397c0fp-12),
		std::make_tuple(-0xc.1eeeb51348828p-7, 0xc.8804209dd17f8p-12),
		std::make_tuple(-0xc.5a2bd2e1bbep-7, 0xc.733d73af55638p-12),
		std::make_tuple(-0xc.4ec3d56605828p-7, 0xc.3407b1f7e5e8p-12),
		std::make_tuple(-0xc.5f9cf4f699ce8p-7, 0xb.fd663100706a8p-12),
		std::make_tuple(-0xa.f6716a67dba38p-7, 0xb.07299ec0060ap-12),
		std::make_tuple(-0x9.d1b8d9cc26cep-7, 0xa.29e5ad37d337p-12),
		std::make_tuple(-0x9.e55e408d6542p-7, 0x9.e05dafff552a8p-12),
		std::make_tuple(-0x8.519fd06204ddp-7, 0x8.b9650f0f51acp-12),
		std::make_tuple(-0xf.cb63af4f243d8p-8, 0x8.210937bb4c6ep-12),
		std::make_tuple(-0xe.5b591b96c887p-8, 0xe.ac6f5d7316b28p-13),
		std::make_tuple(-0xc.9bde6affd0268p-8, 0xc.d9dbff1601be8p-13),
		std::make_tuple(-0xb.2b169c1ffd83p-8, 0xb.148d79cef45c8p-13),
		std::make_tuple(-0x8.80f2a8f5688p-8, 0x8.a50990d6faa8p-13),
		std::make_tuple(-0xc.646a39c9887dp-9, 0xc.a7d6c5340d35p-14),
		std::make_tuple(-0xe.19eaf81fd22e8p-10, 0xf.07c95978714ep-15),
		std::make_tuple(-0xa.7f7a2c2b7e538p-11, 0xc.0d23321d61828p-16),
		std::make_tuple( 0xa.b9750b0852f1p-11, 0xc.351237c4a1288p-16),
		std::make_tuple( 0xf.164dca37b4498p-10, 0xf.85733b90d38a8p-15),
		std::make_tuple( 0xc.84760644e3ea8p-9, 0xc.b308100c2d248p-14),
		std::make_tuple( 0x8.b2bbdbcbba77p-8, 0x8.bc7ef383168b8p-13),
		std::make_tuple( 0xb.0b843ab9167c8p-8, 0xb.07568cfe5a73p-13),
		std::make_tuple( 0xc.47a0f5b8d2f78p-8, 0xc.ac9a4c0d92398p-13),
		std::make_tuple( 0xe.bb561fd5d0068p-8, 0xe.d971942d077cp-13),
		std::make_tuple( 0x8.3ad478182bd28p-7, 0x8.4b395b60088a8p-12),
		std::make_tuple( 0x9.114cdc65ee9a8p-7, 0x9.1add5e10ac528p-12),
		std::make_tuple( 0x9.f958eaa7b0c4p-7, 0x9.ea3fa796f2938p-12),
		std::make_tuple( 0xa.4bdde5165e308p-7, 0xa.666ae06b3bfd8p-12),
		std::make_tuple( 0xb.6306ab51bba08p-7, 0xb.3dc1838312cf8p-12),
		std::make_tuple( 0xb.b88487fb17f5p-7, 0xb.ac4f5e95244cp-12),
		std::make_tuple( 0xc.87592419970d8p-7, 0xc.507173874b078p-12),
		std::make_tuple( 0xc.7d8d18f980eep-7, 0xc.851be4fabc578p-12),
		std::make_tuple( 0xc.adcd3196ebc4p-7, 0xc.cfcdc550c70a8p-12),
		std::make_tuple( 0xc.e2775f65b34e8p-7, 0xd.198558b2efed8p-12),
		std::make_tuple( 0xd.e010759417c18p-7, 0xd.c1c9fb3419f18p-12),
		std::make_tuple( 0xd.c338ec6a34f8p-7, 0xd.db76c7e889e7p-12),
		std::make_tuple( 0xd.9009290a55d7p-7, 0xd.e582f20872fe8p-12),
		std::make_tuple( 0xe.4872901f127c8p-7, 0xe.6370ed4ba505p-12),
		std::make_tuple( 0xe.39de8b66e3e4p-7, 0xe.7a6866caecbp-12),
		std::make_tuple( 0xe.5e0051d38dep-7, 0xe.a8dd91cf76bbp-12),
		std::make_tuple( 0xe.8d26b5f76289p-7, 0xe.daf0c717abc3p-12),
		std::make_tuple( 0xf.40bbc0c14281p-7, 0xf.4db8f1162a608p-12),
		std::make_tuple( 0xe.dc74c588f3fbp-7, 0xf.31f425b1630f8p-12),
		std::make_tuple( 0xf.4d6ab9de38418p-7, 0xf.806f770439a1p-12),
		std::make_tuple( 0xf.84487211f6ddp-7, 0xf.afcfc41e6bcfp-12),
		std::make_tuple( 0x8.2132e75c0fd8p-6, 0x8.106e56262a4bp-11),
		std::make_tuple( 0x8.438ddf975b46p-6, 0x8.2a526fd207518p-11),
		std::make_tuple( 0x8.08737bfdfd998p-6, 0x8.15859b449fec8p-11),
		std::make_tuple( 0x8.19d5766cc49b8p-6, 0x8.25f64530c61cp-11),
		std::make_tuple( 0x8.20499518d71fp-6, 0x8.30af28f4a6f8p-11),
		std::make_tuple( 0x8.b64ea8a68721p-6, 0x8.815c479ab178p-11),
		std::make_tuple( 0x8.0eaeeac394e88p-6, 0x8.356326f113528p-11),
		std::make_tuple( 0x8.0e464b35c2ae8p-6, 0x8.3b83fdc7b67e8p-11),
		std::make_tuple( 0x8.6cd46c94faa3p-6, 0x8.70beb4198b98p-11),
		std::make_tuple( 0x8.8306b28ba7f58p-6, 0x8.818aff2d776p-11),
		std::make_tuple( 0x8.9651612f7619p-6, 0x8.90760d110f8ap-11),
		std::make_tuple( 0x8.5b38528804fa8p-6, 0x8.782e4c74402f8p-11),
		std::make_tuple( 0x8.8bd25d4130fcp-6, 0x8.954b14a28bf18p-11),
		std::make_tuple( 0x8.a69c44b7f98bp-6, 0x8.a74ee7fdaaecp-11),
		std::make_tuple( 0x8.befe6b5cb3f68p-6, 0x8.b7ebf75c54968p-11),
		std::make_tuple( 0x8.667e1bdce0238p-6, 0x8.8ff9b59677dap-11),
		std::make_tuple( 0xf.e08182c2275ap-7, 0x8.5757bbc05b7b8p-11),
		std::make_tuple( 0x8.66467e677e758p-6, 0x8.97bbdda1e16bp-11),
		std::make_tuple( 0x8.8eed6ad22eee8p-6, 0x8.b000ea3b5fbbp-11),
		std::make_tuple( 0x8.6ffad31a783e8p-6, 0x8.a3f7b01cdac88p-11),
		std::make_tuple( 0x8.b77abc45b8af8p-6, 0x8.cb581f50d18f8p-11),
		std::make_tuple( 0x8.db2888fe8e5cp-6, 0x8.e06898bd76b6p-11)
	};
}

template <>
std::vector<std::tuple<long double, long double>> reference_results() {
	return {
		std::make_tuple(-0x9.591830051bfd7efp-6l, 0x9.1e0fc505ad23e91p-11l),
		std::make_tuple(-0x8.db45b19d91c1cc6p-6l, 0x8.dd1e17cb0aa9f87p-11l),
		std::make_tuple(-0x8.c4cb5a7710d01f9p-6l, 0x8.ce82fd5e686a061p-11l),
		std::make_tuple(-0x8.b93b2e9f26ac683p-6l, 0x8.c52802081484c0bp-11l),
		std::make_tuple(-0x8.832c4ddd456e925p-6l, 0x8.a6587d58a845cccp-11l),
		std::make_tuple(-0x9.06de7a096281bcdp-6l, 0x8.e38f934f3f6016ep-11l),
		std::make_tuple(-0x8.e65ae9841ae6a54p-6l, 0x8.cf923e7a71121d2p-11l),
		std::make_tuple(-0x8.9ded366d7105ca8p-6l, 0x8.a779e8890e89654p-11l),
		std::make_tuple(-0x8.af403caa0e6a76fp-6l, 0x8.ab875af7e0144c6p-11l),
		std::make_tuple(-0x8.dc8dff53f7c0db8p-6l, 0x8.bd2416288b5c96dp-11l),
		std::make_tuple(-0x8.f2f5199298ad208p-6l, 0x8.c2f9cd31a5565bbp-11l),
		std::make_tuple(-0x8.8b67a69f4c8ca55p-6l, 0x8.8b012260a5a4a5bp-11l),
		std::make_tuple(-0x8.2ab1d3bf0ed27e3p-6l, 0x8.553bbcae42aa231p-11l),
		std::make_tuple(-0x8.1985f7077c44a80p-6l, 0x8.4702beccff95fb7p-11l),
		std::make_tuple(-0x8.60e183556c4875ap-6l, 0x8.64c6561014693a4p-11l),
		std::make_tuple(-0xf.a320969236cae9fp-7l, 0x8.1660ca7f2e08d2fp-11l),
		std::make_tuple(-0x8.6cebbe36bceccacp-6l, 0x8.5da41cb7ee1c2d7p-11l),
		std::make_tuple(-0xf.fa1ce4eec1487cep-7l, 0x8.1f2303db1b2f646p-11l),
		std::make_tuple(-0x8.0770561d6f86eaep-6l, 0x8.1ce6d8a0a6360b4p-11l),
		std::make_tuple(-0xf.a22af789693df7ap-7l, 0xf.f39ac767f76a5c2p-12l),
		std::make_tuple(-0xf.c32535fa1e6867fp-7l, 0xf.f384cf8af77fb28p-12l),
		std::make_tuple(-0xf.edc2ef78ce726d5p-7l, 0xf.f73786634669789p-12l),
		std::make_tuple(-0xf.b8eba35bbef4325p-7l, 0xf.ca5b8d8693cbbb7p-12l),
		std::make_tuple(-0xf.884ebeb2366b0a4p-7l, 0xf.9d8e6812e734529p-12l),
		std::make_tuple(-0xf.fb33b8f5de31d98p-7l, 0xf.c09cd0c9d459847p-12l),
		std::make_tuple(-0xf.e67081a28304c9fp-7l, 0xf.9e9ca3b517da487p-12l),
		std::make_tuple(-0xe.38cdf0bd20a99b5p-7l, 0xe.affedea00ec4566p-12l),
		std::make_tuple(-0xf.3d207535a71edd3p-7l, 0xf.1808c08f37c2505p-12l),
		std::make_tuple(-0xe.c075c18160c8455p-7l, 0xe.bded27e11fb8774p-12l),
		std::make_tuple(-0xe.140b435cedc23e5p-7l, 0xe.489f241e1b6db06p-12l),
		std::make_tuple(-0xe.9d68ad090417986p-7l, 0xe.6b717b25e7d4edbp-12l),
		std::make_tuple(-0xd.87678ee3e24cf4cp-7l, 0xd.bd8edc23bde28d6p-12l),
		std::make_tuple(-0xd.837eeaccafc1879p-7l, 0xd.944cd80a5a640a8p-12l),
		std::make_tuple(-0xc.9110a1478ef694cp-7l, 0xc.efec35397c0ce8ep-12l),
		std::make_tuple(-0xc.1eeeb513488211ap-7l, 0xc.8804209dd17e04cp-12l),
		std::make_tuple(-0xc.5a2bd2e1bbe0513p-7l, 0xc.733d73af556273ep-12l),
		std::make_tuple(-0xc.4ec3d566058277cp-7l, 0xc.3407b1f7e5e854bp-12l),
		std::make_tuple(-0xc.5f9cf4f699ce7c3p-7l, 0xb.fd663100706f63bp-12l),
		std::make_tuple(-0xa.f6716a67dba2e7cp-7l, 0xb.07299ec0060c76dp-12l),
		std::make_tuple(-0x9.d1b8d9cc26cdf0dp-7l, 0xa.29e5ad37d335487p-12l),
		std::make_tuple(-0x9.e55e408d6541fe3p-7l, 0x9.e05dafff552bc8dp-12l),
		std::make_tuple(-0x8.519fd06204dd0abp-7l, 0x8.b9650f0f51ae898p-12l),
		std::make_tuple(-0xf.cb63af4f243d45ep-8l, 0x8.210937bb4c6d402p-12l),
		std::make_tuple(-0xe.5b591b96c88738fp-8l, 0xe.ac6f5d7316b259ep-13l),
		std::make_tuple(-0xc.9bde6affd026ee5p-8l, 0xc.d9dbff1601bf346p-13l),
		std::make_tuple(-0xb.2b169c1ffd83247p-8l, 0xb.148d79cef45d7e4p-13l),
		std::make_tuple(-0x8.80f2a8f5687fda4p-8l, 0x8.a50990d6faa769bp-13l),
		std::make_tuple(-0xc.646a39c9887d6ebp-9l, 0xc.a7d6c5340d34acdp-14l),
		std::make_tuple(-0xe.19eaf81fd22e9efp-10l, 0xf.07c95978714f12ep-15l),
		std::make_tuple(-0xa.7f7a2c2b7e53f1bp-11l, 0xc.0d23321d6181c20p-16l),
		std::make_tuple( 0xa.b9750b0852ef4b8p-11l, 0xc.351237c4a1283b0p-16l),
		std::make_tuple( 0xf.164dca37b449c9fp-10l, 0xf.85733b90d3895d4p-15l),
		std::make_tuple( 0xc.84760644e3eaeddp-9l, 0xc.b308100c2d25aeep-14l),
		std::make_tuple( 0x8.b2bbdbcbba7763ap-8l, 0x8.bc7ef383168d9dcp-13l),
		std::make_tuple( 0xb.0b843ab9167c5c8p-8l, 0xb.07568cfe5a71d8fp-13l),
		std::make_tuple( 0xc.47a0f5b8d2f794ap-8l, 0xc.ac9a4c0d9235ac1p-13l),
		std::make_tuple( 0xe.bb561fd5d006073p-8l, 0xe.d971942d077d5b7p-13l),
		std::make_tuple( 0x8.3ad478182bd2481p-7l, 0x8.4b395b60088a994p-12l),
		std::make_tuple( 0x9.114cdc65ee9a5afp-7l, 0x9.1add5e10ac52c94p-12l),
		std::make_tuple( 0x9.f958eaa7b0c3c43p-7l, 0x9.ea3fa796f295625p-12l),
		std::make_tuple( 0xa.4bdde5165e30a06p-7l, 0xa.666ae06b3bfd262p-12l),
		std::make_tuple( 0xb.6306ab51bba06c0p-7l, 0xb.3dc1838312ce6c4p-12l),
		std::make_tuple( 0xb.b88487fb17f49e6p-7l, 0xb.ac4f5e95244ba1ep-12l),
		std::make_tuple( 0xc.87592419970d520p-7l, 0xc.507173874b08ac8p-12l),
		std::make_tuple( 0xc.7d8d18f980ee32fp-7l, 0xc.851be4fabc5521ep-12l),
		std::make_tuple( 0xc.adcd3196ebc3f9ap-7l, 0xc.cfcdc550c70c825p-12l),
		std::make_tuple( 0xc.e2775f65b34f018p-7l, 0xd.198558b2efef708p-12l),
		std::make_tuple( 0xd.e010759417c2168p-7l, 0xd.c1c9fb3419efc4bp-12l),
		std::make_tuple( 0xd.c338ec6a34f7fafp-7l, 0xd.db76c7e889e9015p-12l),
		std::make_tuple( 0xd.9009290a55d71c6p-7l, 0xd.e582f20872fe236p-12l),
		std::make_tuple( 0xe.4872901f127ccc5p-7l, 0xe.6370ed4ba505d7fp-12l),
		std::make_tuple( 0xe.39de8b66e3e41ecp-7l, 0xe.7a6866caecb1512p-12l),
		std::make_tuple( 0xe.5e0051d38de02b8p-7l, 0xe.a8dd91cf76b9725p-12l),
		std::make_tuple( 0xe.8d26b5f76288ee1p-7l, 0xe.daf0c717abc353fp-12l),
		std::make_tuple( 0xf.40bbc0c142808d9p-7l, 0xf.4db8f1162a5e2d6p-12l),
		std::make_tuple( 0xe.dc74c588f3fb133p-7l, 0xf.31f425b16311b77p-12l),
		std::make_tuple( 0xf.4d6ab9de3841e20p-7l, 0xf.806f770439a1b24p-12l),
		std::make_tuple( 0xf.84487211f6dd8d6p-7l, 0xf.afcfc41e6bd1622p-12l),
		std::make_tuple( 0x8.2132e75c0fd7f3dp-6l, 0x8.106e56262a4ae55p-11l),
		std::make_tuple( 0x8.438ddf975b46080p-6l, 0x8.2a526fd20753092p-11l),
		std::make_tuple( 0x8.08737bfdfd99117p-6l, 0x8.15859b449fec4f3p-11l),
		std::make_tuple( 0x8.19d5766cc49bab1p-6l, 0x8.25f64530c61d94fp-11l),
		std::make_tuple( 0x8.20499518d71f2b0p-6l, 0x8.30af28f4a6f83d2p-11l),
		std::make_tuple( 0x8.b64ea8a687209a1p-6l, 0x8.815c479ab17856cp-11l),
		std::make_tuple( 0x8.0eaeeac394e8b1cp-6l, 0x8.356326f11351fe4p-11l),
		std::make_tuple( 0x8.0e464b35c2aea6bp-6l, 0x8.3b83fdc7b67cc04p-11l),
		std::make_tuple( 0x8.6cd46c94faa31dap-6l, 0x8.70beb4198b98527p-11l),
		std::make_tuple( 0x8.8306b28ba7f57aep-6l, 0x8.818aff2d776061fp-11l),
		std::make_tuple( 0x8.9651612f7618e6ap-6l, 0x8.90760d110f8a139p-11l),
		std::make_tuple( 0x8.5b38528804fa520p-6l, 0x8.782e4c74402e4aap-11l),
		std::make_tuple( 0x8.8bd25d4130fbf3ap-6l, 0x8.954b14a28bf2180p-11l),
		std::make_tuple( 0x8.a69c44b7f98ad83p-6l, 0x8.a74ee7fdaaeab58p-11l),
		std::make_tuple( 0x8.befe6b5cb3f6a52p-6l, 0x8.b7ebf75c54961f1p-11l),
		std::make_tuple( 0x8.667e1bdce0230e9p-6l, 0x8.8ff9b59677db53bp-11l),
		std::make_tuple( 0xf.e08182c2275aca4p-7l, 0x8.5757bbc05b7c69fp-11l),
		std::make_tuple( 0x8.66467e677e759d3p-6l, 0x8.97bbdda1e16a236p-11l),
		std::make_tuple( 0x8.8eed6ad22eee790p-6l, 0x8.b000ea3b5fbbaa7p-11l),
		std::make_tuple( 0x8.6ffad31a783e261p-6l, 0x8.a3f7b01cdac8cfbp-11l),
		std::make_tuple( 0x8.b77abc45b8af9d2p-6l, 0x8.cb581f50d190080p-11l),
		std::make_tuple( 0x8.db2888fe8e5bd41p-6l, 0x8.e06898bd76b7866p-11l)
	};
}

typedef testing::Types<float, double, long double> MyT;
template <typename T> class DistributionResults : public testing::Test { };
TYPED_TEST_CASE(DistributionResults, MyT);

TYPED_TEST(DistributionResults, CheckPlainIntegration)
{
	using T = TypeParam;

	std::size_t const calls = 100000;

#ifndef HEP_USE_MPI
	auto const result = hep::plain(
#else
	auto const result = hep::mpi_plain(
		MPI_COMM_WORLD,
#endif
		hep::make_integrand<T>(
			arctan<T>,
			1,
			hep::make_dist_params<T>(100, T(-5.0), T(5.0))
		),
		calls
	);

	auto const distribution = result.distributions().at(0).results();
	auto const reference = reference_results<T>();

//	for (auto const dist : distribution)
//	{
//		std::printf(
//			"std::make_tuple(%La, %La),\n",
//			static_cast <long double> (dist.value()),
//			static_cast <long double> (dist.error())
//		);
//	}

	// TODO: check `mid_points()`

	for (std::size_t i = 0; i != 100; ++i)
	{
		EXPECT_DOUBLE_EQ( distribution.at(i).value() , std::get<0>(reference.at(i)) );
		EXPECT_DOUBLE_EQ( distribution.at(i).error() , std::get<1>(reference.at(i)) );
		EXPECT_EQ( distribution.at(i).calls() , calls );
	}
}
