#ifndef SERIALIZER_HPP
#define SERIALIZER_HPP

#include <map>

namespace boost {
    namespace serialization {

        template<template<class T> class SPT>
        class shared_ptr_helper {
            typedef std::map<
                    const void *, // address of object
                    SPT<void> // address shared ptr to single instance
            > object_shared_pointer_map;

            // list of shared_pointers create accessable by raw pointer. This
            // is used to "match up" shared pointers loaded at different
            // points in the archive. Note, we delay construction until
            // it is actually used since this is by default included as
            // a "mix-in" even if shared_ptr isn't used.
            object_shared_pointer_map *m_o_sp;

            struct null_deleter {
                void operator()(void const *) const {}
            };

            template<class Archive, class U>
            friend void boost::serialization::load(
                    Archive &ar,
                    SPT<U> &t,
                    const unsigned int file_version
            );

            struct non_polymorphic {
                template<class U>
                static const boost::serialization::extended_type_info *
                get_object_type(U &) {
                    return &boost::serialization::singleton<
                            typename
                            boost::serialization::type_info_implementation<U>::type
                    >::get_const_instance();
                }
            };

            struct polymorphic {
                template<class U>
                static const boost::serialization::extended_type_info *
                get_object_type(U &u) {
                    return boost::serialization::singleton<
                            typename
                            boost::serialization::type_info_implementation<U>::type
                    >::get_const_instance().get_derived_extended_type_info(u);
                }
            };

        public:
            template<class T>
            void reset(SPT<T> &s, T *t) {
                if (NULL == t) {
                    s.reset();
                    return;
                }
                const boost::serialization::extended_type_info *this_type
                        = &boost::serialization::type_info_implementation<T>::type
                        ::get_const_instance();

                // get pointer to the most derived object's eti.  This is effectively
                // the object type identifer
                typedef typename mpl::if_<
                is_polymorphic < T > ,
                        polymorphic,
                        non_polymorphic
                        > ::type
                type;

                const boost::serialization::extended_type_info *true_type
                        = type::get_object_type(*t);

                // note:if this exception is thrown, be sure that derived pointern
                // is either registered or exported.
                if (NULL == true_type)
                    boost::serialization::throw_exception(
                            boost::archive::archive_exception(
                                    boost::archive::archive_exception::unregistered_class,
                                    this_type->get_debug_info()
                            )
                    );
                // get void pointer to the most derived type
                // this uniquely identifies the object referred to
                // oid = "object identifier"
                const void *oid = void_downcast(
                        *true_type,
                        *this_type,
                        t
                );
                if (NULL == oid)
                    boost::serialization::throw_exception(
                            boost::archive::archive_exception(
                                    boost::archive::archive_exception::unregistered_cast,
                                    true_type->get_debug_info(),
                                    this_type->get_debug_info()
                            )
                    );

                // make tracking array if necessary
                if (NULL == m_o_sp)
                    m_o_sp = new object_shared_pointer_map;

                typename object_shared_pointer_map::iterator i = m_o_sp->find(oid);

                // if it's a new object
                if (i == m_o_sp->end()) {
                    s.reset(t);
                    std::pair<typename object_shared_pointer_map::iterator, bool> result;
                    result = m_o_sp->insert(std::make_pair(oid, s));
                    BOOST_ASSERT(result.second);
                }
                    // if the object has already been seen
                else {
                    s = SPT<T>(i->second, t);
                }
            }

            shared_ptr_helper() :
                    m_o_sp(NULL) {}

            virtual ~shared_ptr_helper() {
                if (NULL != m_o_sp)
                    delete m_o_sp;
            }
        };

    } // namespace serialization
} // namespace boost

namespace boost {
    namespace serialization {

        struct null_deleter {
            void operator()(void const *) const {}
        };

        void *const shared_ptr_helper_id = 0;

        template<class Archive, class T>
        inline void save(
                Archive &ar,
                const std::shared_ptr <T> &t,
                const unsigned int /* file_version */
        ) {
            const T *t_ptr = t.get();
            ar << boost::serialization::make_nvp("px", t_ptr);
        }


        template<class Archive, class T>
        inline void load(
                Archive &ar,
                std::shared_ptr <T> &t,
                const unsigned int /*file_version*/
        ) {
            // The most common cause of trapping here would be serializing
            // something like shared_ptr<int>.  This occurs because int
            // is never tracked by default.  Wrap int in a trackable type
            BOOST_STATIC_ASSERT((tracking_level<T>::value != track_never));
            T *r;
            ar >> boost::serialization::make_nvp("px", r);

            boost::serialization::shared_ptr_helper<std::shared_ptr> &h =
                    ar.template get_helper<shared_ptr_helper<std::shared_ptr> >(
                            shared_ptr_helper_id
                    );
            h.reset(t, r);
        }

        template<class Archive, class T>
        inline void serialize(
                Archive &ar,
                std::shared_ptr <T> &t,
                const unsigned int file_version
        ) {
            // correct shared_ptr serialization depends upon object tracking
            // being used.
            BOOST_STATIC_ASSERT(
                    boost::serialization::tracking_level<T>::value
                    != boost::serialization::track_never
            );
            boost::serialization::split_free(ar, t, file_version);
        }

    }
}
#endif
